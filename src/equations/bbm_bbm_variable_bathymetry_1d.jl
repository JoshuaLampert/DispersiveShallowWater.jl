@doc raw"""
    BBMBBMVariableEquations1D(gravity, D)

BBM-BBM (Benjamin–Bona–Mahony) system in one spatial dimension with spatially varying bathymetry. The equations are given by
```math
\begin{aligned}
  \frac{\partial\eta}{\partial t} + \frac{\partial}{\partial x}((\eta + D)v) - \frac{1}{6}\frac{\partial}{\partial x}(D^2\frac{\partial}{\partial t}\frac{\partial}{\partial x}\eta) &= 0,\\
  \frac{\partial v}{\partial t} + g\frac{\partial\eta}{\partial x} + \frac{\partial}{\partial x}\left(\frac{1}{2}v^2\right) - \frac{1}{6}\frac{\partial^2}{\partial x^2}(D^2\frac{\partial}{\partial t}v) &= 0.
\end{aligned}
```
The unknown quantities of the BBM-BBM equations are the total water height ``\eta`` and the velocity ``v``.
The gravitational constant is denoted by `g` and the bottom topography (bathymetry) ``b = -D > 0``. The water height above the bathymetry is therefore given by
``h = \eta + D``.

One reference for the BBM-BBM system with spatially variying bathymetry can be found in
- Samer Israwi, Henrik Kalisch, Theodoros Katsaounis, Dimitrios Mitsotakis (2022)
  A regularized shallow-water waves system with slip-wall boundary conditions in a basin: theory and numerical analysis
  [DOI: 10.1088/1361-6544/ac3c29](https://doi.org/10.1088/1361-6544/ac3c29)

"""
struct BBMBBMVariableEquations1D{RealT <: Real} <: AbstractBBMBBMEquations{1, 3}
  gravity::RealT # gravitational constant
end

function BBMBBMVariableEquations1D(; gravity_constant)
  BBMBBMVariableEquations1D(gravity_constant)
end

varnames(::BBMBBMVariableEquations1D) = ("eta", "v", "D")

# TODO: Initial condition should not get a `mesh`
"""
    initial_condition_convergence_test(x, t, equations::BBMBBMVariableEquations1D, mesh)

A travelling-wave solution used for convergence tests in a periodic domain.
The bathymetry is constant.

For details see Example 5 in Section 3 from (here adapted for dimensional equations):
- Min Chen (1997)
  Exact Traveling-Wave Solutions to Bidirectional Wave Equations
  [DOI: 10.1023/A:1026667903256](https://doi.org/10.1023/A:1026667903256)
"""
function initial_condition_convergence_test(x,
                                            t,
                                            equations::BBMBBMVariableEquations1D,
                                            mesh)
  g = equations.gravity
  D = 2.0 # constant bathymetry in this case
  c = 5 / 2
  rho = 18 / 5 * sqrt(D * g)
  x_t = mod(x - c * t - xmin(mesh), xmax(mesh) - xmin(mesh)) + xmin(mesh)

  b = 0.5 * sqrt(rho) * x_t / D
  eta = -D + c^2 * rho^2 / (81 * g) +
        5 * c^2 * rho^2 / (108 * g) * (2 / cosh(b)^2 - 3 / cosh(b)^4)
  v = c * (1 - 5 * rho / 18) + 5 * c * rho / 6 / cosh(b)^2
  return SVector(eta, v, D)
end

# TODO: Initial condition should not get a `mesh`
"""
    initial_condition_sin_bathymetry(x, t, equations::BBMBBMVariableEquations1D, mesh)

An initial condition with a gaussion bump as initial water height with still water and
a sine-shaped bathymetry.
"""
function initial_condition_sin_bathymetry(x, t, equations::BBMBBMVariableEquations1D, mesh)
  eta = 1.0 + 2.0 * exp(-12.0 * x^2)
  v = 0.0
  D = -1.0 + 0.1 * sinpi(2.0 * x)
  return SVector(eta, v, D)
end

function create_cache(mesh,
                      equations::BBMBBMVariableEquations1D,
                      solver::Solver,
                      initial_condition,
                      RealT,
                      uEltype)
  #  Assume D is independent of time and compute D evaluated at mesh points once.
  D = Array{RealT}(undef, nnodes(mesh))
  x = grid(solver)
  for i in eachnode(solver)
    D[i] = initial_condition(x[i], 0.0, equations, mesh)[3]
  end
  K = spdiagm(0 => D .^ 2)
  invImDKD_D = (I - 1 / 6 * sparse(solver.D1) * K * sparse(solver.D1)) \ Matrix(solver.D1)
  invImD2K_D = (I - 1 / 6 * sparse(solver.D2) * K) \ Matrix(solver.D1)
  tmp1 = Array{RealT}(undef, nnodes(mesh))
  return (invImDKD_D = invImDKD_D, invImD2K_D = invImD2K_D, tmp1 = tmp1)
end

function create_cache(mesh,
                      equations::BBMBBMVariableEquations1D,
                      solver::UpwindSolver,
                      initial_condition,
                      RealT,
                      uEltype)
  #  Assume D is independent of time and compute D evaluated at mesh points once.
  D = Array{RealT}(undef, nnodes(mesh))
  x = grid(solver)
  for i in eachnode(solver)
    D[i] = initial_condition(x[i], 0.0, equations, mesh)[3]
  end
  K = spdiagm(0 => D .^ 2)
  invImDKD_D = (I - 1 / 6 * sparse(solver.D_min) * K * sparse(solver.D_pl)) \
               Matrix(solver.D_pl)
  invImD2K_D = (I - 1 / 6 * sparse(solver.D2) * K) \ Matrix(solver.D_pl)
  tmp1 = Array{RealT}(undef, nnodes(mesh))
  return (invImDKD_D = invImDKD_D, invImD2K_D = invImD2K_D, tmp1 = tmp1)
end

# Discretization that conserves the mass (for eta and u) and the energy for periodic boundary conditions, see
# - Hendrik Ranocha, Dimitrios Mitsotakis and David I. Ketcheson (2020)
#   A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations
#   [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)
# Here, adapted for spatially varying bathymetry.
function rhs!(du_ode, u_ode, t, mesh, equations::BBMBBMVariableEquations1D,
              initial_condition,
              ::BoundaryConditionPeriodic, solver, cache)
  @unpack invImDKD_D, invImD2K_D, tmp1 = cache

  u = wrap_array(u_ode, mesh, equations, solver)
  du = wrap_array(du_ode, mesh, equations, solver)

  eta = view(u, 1, :)
  v = view(u, 2, :)
  D = view(u, 3, :)
  deta = view(du, 1, :)
  dv = view(du, 2, :)
  dD = view(du, 3, :)
  fill!(dD, zero(eltype(dD)))

  @. tmp1 = -(D * v + eta * v)
  mul!(deta, invImDKD_D, tmp1)

  @. tmp1 = -(equations.gravity * eta + 0.5 * v^2)
  mul!(dv, invImD2K_D, tmp1)

  return nothing
end

@inline function waterheight_total(u, equations::BBMBBMVariableEquations1D)
  return u[1]
end

@inline function velocity(u, equations::BBMBBMVariableEquations1D)
  return u[2]
end

@inline function bathymetry(u, equations::BBMBBMVariableEquations1D)
  return u[3]
end

@inline function waterheight(u, equations::BBMBBMVariableEquations1D)
  return waterheight_total(u, equations) + bathymetry(u, equations)
end

@inline function energy_total(u, equations::BBMBBMVariableEquations1D)
  eta, v, D = u
  e = equations.gravity * eta^2 + (D + eta) * v^2
  return e
end

@inline entropy(u, equations::BBMBBMVariableEquations1D) = energy_total(u, equations)
