@doc raw"""
    BBMBBMEquations1D(gravity, D)

BBM-BBM (Benjamin–Bona–Mahony) system in one spatial dimension. The equations are given by
```math
\begin{aligned}
  \frac{\partial\eta}{\partial t} + \frac{\partial}{\partial x}((\eta + D)v) - \frac{1}{6}D^2\frac{\partial}{\partial t}\frac{\partial^2}{\partial x^2}\eta &= 0,\\
  \frac{\partial v}{\partial t} + g\frac{\partial\eta}{\partial x} + \frac{\partial}{\partial x}\left(\frac{1}{2}v^2\right) - \frac{1}{6}D^2\frac{\partial}{\partial t}\frac{\partial^2}{\partial x^2}v &= 0.
\end{aligned}
```
The unknown quantities of the BBM-BBM equations are the total water height ``\eta`` and the velocity ``v``.
The gravitational constant is denoted by `g` and the constant bottom topography (bathymetry) ``b = -D``. The water height above the bathymetry is therefore given by
``h = \eta + D``.

One reference for the BBM-BBM system can be found in
- Jerry L. Bona, Min Chen (1998)
  A Boussinesq system for two-way propagation of nonlinear dispersive waves
  [DOI: 10.1016/S0167-2789(97)00249-2](https://doi.org/10.1016/S0167-2789(97)00249-2)

"""
struct BBMBBMEquations1D{RealT <: Real} <: AbstractBBMBBMEquations{1, 2}
    gravity::RealT # gravitational constant
    D::RealT       # constant bathymetry
end

function BBMBBMEquations1D(; gravity_constant, D = 1.0)
    BBMBBMEquations1D(gravity_constant, D)
end

varnames(::BBMBBMEquations1D) = ("eta", "v")

# TODO: Initial condition should not get a `mesh`
"""
    initial_condition_convergence_test(x, t, equations::BBMBBMEquations1D, mesh)

A travelling-wave solution used for convergence tests in a periodic domain.

For details see Example 5 in Section 3 from (here adapted for dimensional equations):
- Min Chen (1997)
  Exact Traveling-Wave Solutions to Bidirectional Wave Equations
  [DOI: 10.1023/A:1026667903256](https://doi.org/10.1023/A:1026667903256)
"""
function initial_condition_convergence_test(x, t, equations::BBMBBMEquations1D, mesh)
    g = equations.gravity
    c = 5 / 2 * sqrt(equations.D * g)
    rho = 18 / 5
    x_t = mod(x - c * t - xmin(mesh), xmax(mesh) - xmin(mesh)) + xmin(mesh)

    theta = 0.5 * sqrt(rho) * x_t / equations.D
    eta = -equations.D + c^2 * rho^2 / (81 * g) +
          5 * c^2 * rho^2 / (108 * g) * (2 * sech(theta)^2 - 3 * sech(theta)^4)
    v = c * (1 - 5 * rho / 18) + 5 * c * rho / 6 * sech(theta)^2
    return SVector(eta, v)
end

function create_cache(mesh,
                      equations::BBMBBMEquations1D,
                      solver,
                      initial_condition,
                      RealT,
                      uEltype)
    invImD2_D = (I - 1 / 6 * equations.D^2 * sparse(solver.D2)) \ Matrix(solver.D1)
    tmp1 = Array{RealT}(undef, nnodes(mesh))
    return (invImD2_D = invImD2_D, tmp1 = tmp1)
end

# Discretization that conserves the mass (for eta and v) and the energy for periodic boundary conditions, see
# - Hendrik Ranocha, Dimitrios Mitsotakis and David I. Ketcheson (2020)
#   A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations
#   [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)
function rhs!(du_ode, u_ode, t, mesh, equations::BBMBBMEquations1D, initial_condition,
              ::BoundaryConditionPeriodic, solver, cache)
    @unpack invImD2_D, tmp1 = cache

    q = wrap_array(u_ode, mesh, equations, solver)
    dq = wrap_array(du_ode, mesh, equations, solver)

    eta = view(q, 1, :)
    v = view(q, 2, :)
    deta = view(dq, 1, :)
    dv = view(dq, 2, :)

    # energy and mass conservative semidiscretization
    @. tmp1 = -(equations.D * v + eta * v)
    mul!(deta, invImD2_D, tmp1)

    @. tmp1 = -(equations.gravity * eta + 0.5 * v^2)
    mul!(dv, invImD2_D, tmp1)

    return nothing
end

@inline function waterheight_total(q, equations::BBMBBMEquations1D)
    return q[1]
end

@inline function velocity(q, equations::BBMBBMEquations1D)
    return q[2]
end

@inline function bathymetry(q, equations::BBMBBMEquations1D)
    return -equations.D
end

@inline function waterheight(q, equations::BBMBBMEquations1D)
    return waterheight_total(q, equations) - bathymetry(q, equations)
end

@inline function energy_total(q, equations::BBMBBMEquations1D)
    eta, v = q
    e = 0.5 * (equations.gravity * eta^2 + (equations.D + eta) * v^2)
    return e
end

@inline entropy(q, equations::BBMBBMEquations1D) = energy_total(q, equations)
