@doc raw"""
    BBMBBMEquations1D(gravity, D)

BBM-BBM (Benjamin–Bona–Mahony) system in one spatial dimension. The equations are given by
```math
\begin{aligned}
  \eta_t + ((\eta + D)v)_x - \frac{1}{6}D^2\eta_{xxt} &= 0,\\
  v_t + g\eta_x + \left(\frac{1}{2}v^2\right)_x - \frac{1}{6}D^2v_{xxt} &= 0.
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

varnames(::typeof(prim2prim), ::BBMBBMEquations1D) = ("η", "v")
varnames(::typeof(prim2cons), ::BBMBBMEquations1D) = ("h", "hv")

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

"""
    initial_condition_manufactured(x, t, equations::BBMBBMEquations1D, mesh)

A smooth manufactured solution in combination with [`source_terms_manufactured`](@ref).
"""
function initial_condition_manufactured(x, t,
                                        equations::BBMBBMEquations1D,
                                        mesh)
    eta = exp(t) * cospi(2 * (x - 2 * t))
    v = exp(t / 2) * sinpi(2 * (x - t / 2))
    return SVector(eta, v)
end

"""
    source_terms_manufactured(q, x, t, equations::BBMBBMEquations1D, mesh)

A smooth manufactured solution in combination with [`initial_condition_manufactured`](@ref).
"""
function source_terms_manufactured(q, x, t, equations::BBMBBMEquations1D)
    g = equations.gravity
    D = equations.D
    a3 = cospi(t - 2 * x)
    a4 = sinpi(t - 2 * x)
    a5 = sinpi(2 * t - 4 * x)
    a6 = sinpi(4 * t - 2 * x)
    a7 = cospi(4 * t - 2 * x)
    dq1 = -2 * pi^2 * D^2 * (4 * pi * a6 - a7) * exp(t) / 3 + 2 * pi * D * exp(t / 2) * a3 -
          2 * pi * exp(3 * t / 2) * a4 * a6 + 2 * pi * exp(3 * t / 2) * a3 * a7 -
          4 * pi * exp(t) * a6 + exp(t) * a7
    dq2 = -pi^2 * D^2 * (a4 + 2 * pi * a3) * exp(t / 2) / 3 + 2 * pi * g * exp(t) * a6 -
          exp(t / 2) * a4 / 2 - pi * exp(t / 2) * a3 - pi * exp(t) * a5

    return SVector(dq1, dq2)
end
#
function create_cache(mesh,
                      equations::BBMBBMEquations1D,
                      solver,
                      initial_condition,
                      RealT,
                      uEltype)
    if solver.D1 isa PeriodicDerivativeOperator ||
       solver.D1 isa UniformPeriodicCoupledOperator
        invImD2 = inv(I - 1 / 6 * equations.D^2 * Matrix(solver.D2))
    elseif solver.D1 isa PeriodicUpwindOperators
        invImD2 = inv(I - 1 / 6 * equations.D^2 * Matrix(solver.D2))
    else
        @error "unknown type of first-derivative operator"
    end
    tmp1 = Array{RealT}(undef, nnodes(mesh)) # tmp1 is needed for the `RelaxationCallback`
    return (invImD2 = invImD2, tmp1 = tmp1)
end

# Discretization that conserves the mass (for eta and v) and the energy for periodic boundary conditions, see
# - Hendrik Ranocha, Dimitrios Mitsotakis and David I. Ketcheson (2020)
#   A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations
#   [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)
function rhs!(du_ode, u_ode, t, mesh, equations::BBMBBMEquations1D, initial_condition,
              ::BoundaryConditionPeriodic, source_terms, solver, cache)
    @unpack invImD2 = cache

    q = wrap_array(u_ode, mesh, equations, solver)
    dq = wrap_array(du_ode, mesh, equations, solver)

    eta = view(q, 1, :)
    v = view(q, 2, :)
    deta = view(dq, 1, :)
    dv = view(dq, 2, :)

    # energy and mass conservative semidiscretization
    if solver.D1 isa PeriodicDerivativeOperator ||
       solver.D1 isa UniformPeriodicCoupledOperator
        deta[:] = -solver.D1 * (equations.D * v + eta .* v)
        dv[:] = -solver.D1 * (equations.gravity * eta + 0.5 * v .^ 2)
    elseif solver.D1 isa PeriodicUpwindOperators
        deta[:] = -solver.D1.central * (equations.D * v + eta .* v)
        dv[:] = -solver.D1.central * (equations.gravity * eta + 0.5 * v .^ 2)
    else
        @error "unknown type of first-derivative operator"
    end

    calc_sources!(dq, q, t, source_terms, equations, solver)

    deta[:] = invImD2 * deta
    dv[:] = invImD2 * dv

    return nothing
end

@inline function prim2cons(q, equations::BBMBBMEquations1D)
    eta, v = q

    h = eta + equations.D
    hv = h * v
    return SVector(h, hv)
end

@inline function cons2prim(u, equations::BBMBBMEquations1D)
    h, hv = u

    eta = h - equations.D
    v = hv / h
    return SVector(eta, v)
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
