@doc raw"""
    BBMBBMEquations1D(gravity, D, eta0 = 0.0)

BBM-BBM (Benjamin–Bona–Mahony) system in one spatial dimension. The equations are given by
```math
\begin{aligned}
  \eta_t + ((\eta + D)v)_x - \frac{1}{6}D^2\eta_{xxt} &= 0,\\
  v_t + g\eta_x + \left(\frac{1}{2}v^2\right)_x - \frac{1}{6}D^2v_{xxt} &= 0.
\end{aligned}
```
The unknown quantities of the BBM-BBM equations are the total water height ``\eta`` and the velocity ``v``.
The gravitational constant is denoted by `g` and the constant bottom topography (bathymetry) ``b = \eta_0 - D``. The water height above the bathymetry is therefore given by
``h = \eta - \eta_0 + D``. The BBM-BBM equations are only implemented for ``\eta_0 = 0``.

One reference for the BBM-BBM system can be found in
- Jerry L. Bona, Min Chen (1998)
  A Boussinesq system for two-way propagation of nonlinear dispersive waves
  [DOI: 10.1016/S0167-2789(97)00249-2](https://doi.org/10.1016/S0167-2789(97)00249-2)

"""
struct BBMBBMEquations1D{RealT <: Real} <: AbstractBBMBBMEquations{1, 2}
    gravity::RealT # gravitational constant
    D::RealT       # constant bathymetry
    eta0::RealT    # constant still-water surface
end

function BBMBBMEquations1D(; gravity_constant, D = 1.0, eta0 = 0.0)
    eta0 == 0.0 || @warn "The still-water surface needs to be 0 for the BBM-BBM equations"
    BBMBBMEquations1D(gravity_constant, D, eta0)
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

"""
    initial_condition_manufactured_reflecting(x, t, equations::BBMBBMEquations1D, mesh)

A smooth manufactured solution for reflecting boundary conditions in combination
with [`source_terms_manufactured_reflecting`](@ref).
"""
function initial_condition_manufactured_reflecting(x, t,
                                                   equations::BBMBBMEquations1D,
                                                   mesh)
    eta = exp(2 * t) * cospi(x)
    v = exp(t) * x * sinpi(x)
    return SVector(eta, v)
end

"""
    source_terms_manufactured_reflecting(q, x, t, equations::BBMBBMEquations1D, mesh)

A smooth manufactured solution for reflecting boundary conditions in combination
with [`initial_condition_manufactured_reflecting`](@ref).
"""
function source_terms_manufactured_reflecting(q, x, t, equations::BBMBBMEquations1D)
    g = equations.gravity
    D = equations.D
    a1 = cospi(2 * x)
    a2 = sinpi(2 * x)
    a8 = cospi(x)
    a9 = sinpi(x)
    a10 = exp(t)
    a11 = exp(2 * t)
    dq1 = (pi^2 * D^2 * a10 * a8 / 3 + pi * D * x * a8 + D * a9 + pi * x * a11 * a1 +
           a11 * a2 / 2 + 2 * a10 * a8) * a10
    dq2 = (pi * D^2 * (pi * x * a9 - 2 * a8) / 6 - pi * g * a10 * a9 +
           pi * x^2 * a10 * a2 / 2 + x * a10 * a9^2 + x * a9) * a10

    return SVector(dq1, dq2)
end

function create_cache(mesh, equations::BBMBBMEquations1D,
                      solver, initial_condition,
                      ::BoundaryConditionPeriodic,
                      RealT, uEltype)
    D = equations.D
    invImD2 = lu(I - 1 / 6 * D^2 * sparse(solver.D2))
    tmp2 = Array{RealT}(undef, nnodes(mesh))
    tmp3 = similar(tmp2)
    tmp4 = similar(tmp2)
    return (invImD2 = invImD2, tmp2 = tmp2, tmp3 = tmp3, tmp4 = tmp4)
end

function create_cache(mesh, equations::BBMBBMEquations1D,
                      solver, initial_condition,
                      ::BoundaryConditionReflecting,
                      RealT, uEltype)
    D = equations.D
    N = nnodes(mesh)
    M = mass_matrix(solver.D1)
    Pd = BandedMatrix((-1 => fill(one(real(mesh)), N - 2),), (N, N - 2))
    D2d = (sparse(solver.D2) * Pd)[2:(end - 1), :]
    # homogeneous Dirichlet boundary conditions
    invImD2d = lu(I - 1 / 6 * D^2 * D2d)
    m = diag(M)
    m[1] = 0
    m[end] = 0
    PdM = Diagonal(m)

    # homogeneous Neumann boundary conditions
    if solver.D1 isa DerivativeOperator ||
       solver.D1 isa UniformCoupledOperator
        D1_b = BandedMatrix(solver.D1)
        invImD2n = lu(I + 1 / 6 * D^2 * inv(M) * D1_b' * PdM * D1_b)
    elseif solver.D1 isa UpwindOperators
        D1plus_b = BandedMatrix(solver.D1.plus)
        invImD2n = lu(I + 1 / 6 * D^2 * inv(M) * D1plus_b' * PdM * D1plus_b)
    else
        @error "unknown type of first-derivative operator: $(typeof(solver.D1))"
    end
    tmp2 = Array{RealT}(undef, nnodes(mesh))
    tmp3 = similar(tmp2)
    tmp4 = Array{RealT}(undef, nnodes(mesh) - 2)
    return (invImD2d = invImD2d, invImD2n = invImD2n, tmp2 = tmp2, tmp3 = tmp3, tmp4 = tmp4)
end

# Discretization that conserves the mass (for eta and v) and the energy for periodic boundary conditions, see
# - Hendrik Ranocha, Dimitrios Mitsotakis and David I. Ketcheson (2020)
#   A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations
#   [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)
function rhs!(du_ode, u_ode, t, mesh, equations::BBMBBMEquations1D, initial_condition,
              ::BoundaryConditionPeriodic, source_terms, solver, cache)
    @unpack invImD2, tmp1, tmp2, tmp3, tmp4 = cache

    q = wrap_array(u_ode, mesh, equations, solver)
    dq = wrap_array(du_ode, mesh, equations, solver)

    eta = view(q, 1, :)
    v = view(q, 2, :)
    deta = view(dq, 1, :)
    dv = view(dq, 2, :)

    D = equations.D
    # energy and mass conservative semidiscretization
    if solver.D1 isa PeriodicDerivativeOperator ||
       solver.D1 isa UniformPeriodicCoupledOperator
        @timeit timer() "deta hyperbolic" deta[:]=-solver.D1 * (D * v + eta .* v)
        @timeit timer() "dv hyperbolic" dv[:]=-solver.D1 *
                                              (equations.gravity * eta + 0.5 * v .^ 2)
    elseif solver.D1 isa PeriodicUpwindOperators
        # Note that the upwind operators here are not actually used
        # We would need to define two different matrices `invImD2` for eta and v for energy conservation
        # To really use the upwind operators, we can use them with `BBMBBMVariableEquations1D`
        @timeit timer() "deta hyperbolic" deta[:]=-solver.D1.central * (D * v + eta .* v)
        @timeit timer() "dv hyperbolic" dv[:]=-solver.D1.central *
                                              (equations.gravity * eta + 0.5 * v .^ 2)
    else
        @error "unknown type of first-derivative operator: $(typeof(solver.D1))"
    end

    @timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations, solver)

    # To use the in-place version `ldiv!` instead of `\`, we need temporary arrays
    # since `deta` and `dv` are not stored contiguously
    @timeit timer() "deta elliptic" begin
        tmp1[:] = deta
        ldiv!(tmp3, invImD2, tmp1)
        deta[:] = tmp3
    end
    @timeit timer() "dv elliptic" begin
        tmp2[:] = dv
        ldiv!(tmp4, invImD2, tmp2)
        dv[:] = tmp4
    end
    return nothing
end

# Discretization that conserves the mass (for eta) and the energy for periodic boundary conditions, see
# - Hendrik Ranocha, Dimitrios Mitsotakis and David I. Ketcheson (2020)
#   A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations
#   [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)
function rhs!(du_ode, u_ode, t, mesh, equations::BBMBBMEquations1D, initial_condition,
              ::BoundaryConditionReflecting, source_terms, solver, cache)
    @unpack invImD2d, invImD2n, tmp1, tmp2, tmp3, tmp4 = cache

    q = wrap_array(u_ode, mesh, equations, solver)
    dq = wrap_array(du_ode, mesh, equations, solver)

    eta = view(q, 1, :)
    v = view(q, 2, :)
    deta = view(dq, 1, :)
    dv = view(dq, 2, :)

    D = equations.D
    # energy and mass conservative semidiscretization
    if solver.D1 isa DerivativeOperator ||
       solver.D1 isa UniformCoupledOperator
        @timeit timer() "deta hyperbolic" deta[:]=-solver.D1 * (D * v + eta .* v)
        @timeit timer() "dv hyperbolic" dv[:]=-solver.D1 *
                                              (equations.gravity * eta + 0.5 * v .^ 2)
    elseif solver.D1 isa UpwindOperators
        @timeit timer() "deta hyperbolic" deta[:]=-solver.D1.minus * (D * v + eta .* v)
        @timeit timer() "dv hyperbolic" dv[:]=-solver.D1.plus *
                                              (equations.gravity * eta + 0.5 * v .^ 2)
    else
        @error "unknown type of first-derivative operator: $(typeof(solver.D1))"
    end

    @timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations, solver)

    # To use the in-place version `ldiv!` instead of `\`, we need temporary arrays
    # since `deta` and `dv` are not stored contiguously
    @timeit timer() "deta elliptic" begin
        tmp1[:] = deta
        ldiv!(tmp3, invImD2n, tmp1)
        deta[:] = tmp3
    end
    @timeit timer() "dv elliptic" begin
        tmp2[:] = dv
        ldiv!(tmp4, invImD2d, tmp2[2:(end - 1)])
        dv[1] = dv[end] = zero(eltype(dv))
        dv[2:(end - 1)] = tmp4
    end

    return nothing
end

@inline function prim2cons(q, equations::BBMBBMEquations1D)
    eta, v = q

    h = eta - equations.eta0 + equations.D
    hv = h * v
    return SVector(h, hv)
end

@inline function cons2prim(u, equations::BBMBBMEquations1D)
    h, hv = u

    eta = h + equations.eta0 - equations.D
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
    return equations.eta0 - equations.D
end

@inline function waterheight(q, equations::BBMBBMEquations1D)
    return waterheight_total(q, equations) - bathymetry(q, equations)
end

@inline function energy_total(q, equations::BBMBBMEquations1D)
    eta, v = q
    D = still_waterdepth(q, equations)
    e = 0.5 * (equations.gravity * eta^2 + (D + eta - equations.eta0) * v^2)
    return e
end

@inline entropy(q, equations::BBMBBMEquations1D) = energy_total(q, equations)
