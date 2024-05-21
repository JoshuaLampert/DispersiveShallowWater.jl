@doc raw"""
    BBMBBMVariableEquations1D(gravity, eta0 = 1.0)

BBM-BBM (Benjamin–Bona–Mahony) system in one spatial dimension with spatially varying bathymetry. The equations are given by
```math
\begin{aligned}
  \eta_t + ((\eta + D)v)_x - \frac{1}{6}(D^2\eta_{xt})_x &= 0,\\
  v_t + g\eta_x + \left(\frac{1}{2}v^2\right)_x - \frac{1}{6}(D^2v_t)_{xx} &= 0.
\end{aligned}
```
The unknown quantities of the BBM-BBM equations are the total water height ``\eta`` and the velocity ``v``.
The gravitational constant is denoted by `g` and the bottom topography (bathymetry) ``b = \eta_0 - D``. The water height above the bathymetry is therefore given by
``h = \eta - \eta_0 + D``. The BBM-BBM equations are only implemented for ``\eta_0 = 0``.

One reference for the BBM-BBM system with spatially varying bathymetry can be found in
- Samer Israwi, Henrik Kalisch, Theodoros Katsaounis, Dimitrios Mitsotakis (2022)
  A regularized shallow-water waves system with slip-wall boundary conditions in a basin: theory and numerical analysis
  [DOI: 10.1088/1361-6544/ac3c29](https://doi.org/10.1088/1361-6544/ac3c29)

"""
struct BBMBBMVariableEquations1D{RealT <: Real} <: AbstractBBMBBMEquations{1, 3}
    gravity::RealT # gravitational constant
    eta0::RealT    # constant still-water surface
end

function BBMBBMVariableEquations1D(; gravity_constant, eta0 = 0.0)
    eta0 == 0.0 || @warn "The still-water surface needs to be 0 for the BBM-BBM equations"
    BBMBBMVariableEquations1D(gravity_constant, eta0)
end

varnames(::typeof(prim2prim), ::BBMBBMVariableEquations1D) = ("η", "v", "D")
varnames(::typeof(prim2cons), ::BBMBBMVariableEquations1D) = ("h", "hv", "b")

"""
    initial_condition_convergence_test(x, t, equations::BBMBBMVariableEquations1D, mesh)

A travelling-wave solution used for convergence tests in a periodic domain.
The bathymetry is constant.

For details see Example 5 in Section 3 from (here adapted for dimensional equations):
- Min Chen (1997)
  Exact Traveling-Wave Solutions to Bidirectional Wave Equations
  [DOI: 10.1023/A:1026667903256](https://doi.org/10.1023/A:1026667903256)
"""
function initial_condition_convergence_test(x, t,
                                            equations::BBMBBMVariableEquations1D,
                                            mesh)
    g = equations.gravity
    D = 2.0 # constant bathymetry in this case
    c = 5 / 2 * sqrt(D * g)
    rho = 18 / 5
    x_t = mod(x - c * t - xmin(mesh), xmax(mesh) - xmin(mesh)) + xmin(mesh)

    theta = 0.5 * sqrt(rho) * x_t / D
    eta = -D + c^2 * rho^2 / (81 * g) +
          5 * c^2 * rho^2 / (108 * g) * (2 * sech(theta)^2 - 3 * sech(theta)^4)
    v = c * (1 - 5 * rho / 18) + 5 * c * rho / 6 * sech(theta)^2
    return SVector(eta, v, D)
end

"""
    initial_condition_manufactured(x, t, equations::BBMBBMVariableEquations1D, mesh)

A smooth manufactured solution in combination with [`source_terms_manufactured`](@ref).
"""
function initial_condition_manufactured(x, t,
                                        equations::BBMBBMVariableEquations1D,
                                        mesh)
    eta = exp(t) * cospi(2 * (x - 2 * t))
    v = exp(t / 2) * sinpi(2 * (x - t / 2))
    D = 5 + 2 * cospi(2 * x)
    return SVector(eta, v, D)
end

"""
    source_terms_manufactured(q, x, t, equations::BBMBBMVariableEquations1D, mesh)

A smooth manufactured solution in combination with [`initial_condition_manufactured`](@ref).
"""
function source_terms_manufactured(q, x, t, equations::BBMBBMVariableEquations1D)
    g = equations.gravity
    a1 = cospi(2 * x)
    a2 = sinpi(2 * x)
    a3 = cospi(t - 2 * x)
    a4 = sinpi(t - 2 * x)
    a5 = sinpi(2 * t - 4 * x)
    a6 = sinpi(4 * t - 2 * x)
    a7 = cospi(4 * t - 2 * x)
    dq1 = -2 * pi^2 * (4 * pi * a6 - a7) * (2 * a1 + 5)^2 * exp(t) / 3 +
          8 * pi^2 * (a6 + 4 * pi * a7) * (2 * a1 + 5) * exp(t) * a2 / 3 +
          2 * pi * (2 * a1 + 5) * exp(t / 2) * a3 - 2 * pi * exp(3 * t / 2) * a4 * a6 +
          2 * pi * exp(3 * t / 2) * a3 * a7 + 4 * pi * exp(t / 2) * a2 * a4 -
          4 * pi * exp(t) * a6 + exp(t) * a7
    dq2 = 2 * pi * g * exp(t) * a6 -
          pi^2 *
          (8 * (2 * pi * a4 - a3) * (2 * a1 + 5) * a2 +
           (a4 + 2 * pi * a3) * (2 * a1 + 5)^2 +
           4 * (a4 + 2 * pi * a3) * (16 * sinpi(x)^4 - 26 * sinpi(x)^2 + 7)) * exp(t / 2) /
          3 - exp(t / 2) * a4 / 2 - pi * exp(t / 2) * a3 - pi * exp(t) * a5

    return SVector(dq1, dq2, zero(dq1))
end

"""
    initial_condition_manufactured_reflecting(x, t, equations::BBMBBMVariableEquations1D, mesh)

A smooth manufactured solution for reflecting boundary conditions in combination
with [`source_terms_manufactured_reflecting`](@ref).
"""
function initial_condition_manufactured_reflecting(x, t,
                                                   equations::BBMBBMVariableEquations1D,
                                                   mesh)
    eta = exp(2 * t) * cospi(x)
    v = exp(t) * x * sinpi(x)
    D = 5 + 2 * cospi(2 * x)
    return SVector(eta, v, D)
end

"""
    source_terms_manufactured_reflecting(q, x, t, equations::BBMBBMVariableEquations1D, mesh)

A smooth manufactured solution for reflecting boundary conditions in combination
with [`initial_condition_manufactured_reflecting`](@ref).
"""
function source_terms_manufactured_reflecting(q, x, t, equations::BBMBBMVariableEquations1D)
    g = equations.gravity
    a1 = cospi(2 * x)
    a2 = sinpi(2 * x)
    a8 = cospi(x)
    a9 = sinpi(x)
    a10 = exp(t)
    a11 = exp(2 * t)
    a12 = cospi(3*x)
    a13 = sinpi(3*x)
    dq1 = (pi*x*a11*a1 + 4*pi*x*a8 + 3*pi*x*a12 + 20*pi^2*(1 - a1)^2*a10*a8/3 + a11*a2/2 + 2*a10*a8 + 7*pi^2*a10*a8/3 + 14*pi^2*a10*a12 + 4*a9 + a13)*a10
    dq2 = (-pi*g*a10*a9 + pi*x^2*a10*a2/2 + x*a10*a9^2 + x*a9 + pi*(400*pi*x*a9^5 - 824*pi*x*a9^3 + 385*pi*x*a9 - 160*a9^4*a8 + 336*a9^2*a8 - 98*a8)/6)*a10

    return SVector(dq1, dq2, zero(dq1))
end

"""
    initial_condition_dingemans(x, t, equations::BBMBBMVariableEquations1D, mesh)

The initial condition that uses the dispersion relation of the Euler equations
to approximate waves generated by a wave maker as it is done by experiments of
Dingemans. The topography is a trapezoidal.

!!! warning "Translation of water height"
    The initial condition for the water height is translated to be around 0, which is
    needed for the simulation because the `BBMBBMVariableEquations1D` are only implemented
    for ``\eta_0 = 0``.

References:
- Magnus Svärd, Henrik Kalisch (2023)
  A novel energy-bounded Boussinesq model and a well-balanced and stable numerical discretization
  [arXiv: 2302.09924](https://arxiv.org/abs/2302.09924)
- Maarten W. Dingemans (1994)
  Comparison of computations with Boussinesq-like models and laboratory measurements
  [link](https://repository.tudelft.nl/islandora/object/uuid:c2091d53-f455-48af-a84b-ac86680455e9/datastream/OBJ/download)
"""
function initial_condition_dingemans(x, t, equations::BBMBBMVariableEquations1D, mesh)
    h0 = 0.8
    A = 0.02
    # omega = 2*pi/(2.02*sqrt(2))
    k = 0.8406220896381442 # precomputed result of find_zero(k -> omega^2 - equations.gravity * k * tanh(k * h0), 1.0) using Roots.jl
    if x < -30.5 * pi / k || x > -8.5 * pi / k
        h = 0.0
    else
        h = A * cos(k * x)
    end
    v = sqrt(equations.gravity / k * tanh(k * h0)) * h / h0
    if 11.01 <= x && x < 23.04
        b = 0.6 * (x - 11.01) / (23.04 - 11.01)
    elseif 23.04 <= x && x < 27.04
        b = 0.6
    elseif 27.04 <= x && x < 33.07
        b = 0.6 * (33.07 - x) / (33.07 - 27.04)
    else
        b = 0.0
    end
    # Here, we compute eta - h0!! To obtain the original eta, h0 = 0.8 needs to be added again!
    # This is because the BBM-BBM equations are only implemented for eta0 = 0
    eta = h
    D = h0 - b
    return SVector(eta, v, D)
end

function create_cache(mesh, equations::BBMBBMVariableEquations1D,
                      solver, initial_condition,
                      ::BoundaryConditionPeriodic,
                      RealT, uEltype)
    #  Assume D is independent of time and compute D evaluated at mesh points once.
    D = Array{RealT}(undef, nnodes(mesh))
    x = grid(solver)
    for i in eachnode(solver)
        D[i] = still_waterdepth(initial_condition(x[i], 0.0, equations, mesh), equations)
    end
    K = Diagonal(D .^ 2)
    if solver.D1 isa PeriodicDerivativeOperator ||
       solver.D1 isa UniformPeriodicCoupledOperator
        invImDKD = inv(I - 1 / 6 * Matrix(solver.D1) * K * Matrix(solver.D1))
    elseif solver.D1 isa PeriodicUpwindOperators
        invImDKD = inv(I - 1 / 6 * Matrix(solver.D1.minus) * K * Matrix(solver.D1.plus))
    else
        @error "unknown type of first-derivative operator: $(typeof(solver.D1))"
    end
    invImD2K = inv(I - 1 / 6 * Matrix(solver.D2) * K)
    return (invImDKD = invImDKD, invImD2K = invImD2K, D = D)
end

function create_cache(mesh, equations::BBMBBMVariableEquations1D,
                      solver, initial_condition,
                      ::BoundaryConditionReflecting,
                      RealT, uEltype)
    #  Assume D is independent of time and compute D evaluated at mesh points once.
    D = Array{RealT}(undef, nnodes(mesh))
    x = grid(solver)
    for i in eachnode(solver)
        D[i] = still_waterdepth(initial_condition(x[i], 0.0, equations, mesh), equations)
    end
    K = Diagonal(D .^ 2)
    K_i = Diagonal(D[2:(end - 1)] .^ 2)
    N = nnodes(mesh)
    M = mass_matrix(solver.D1)
    Pd = BandedMatrix((-1 => fill(one(real(mesh)), N - 2),), (N, N - 2))
    D2d = (sparse(solver.D2) * Pd)[2:(end - 1), :]
    # homogeneous Dirichlet boundary conditions
    invImD2Kd = inv(I - 1 / 6 * D2d * K_i)
    m = diag(M)
    m[1] = 0
    m[end] = 0
    PdM = Diagonal(m)

    # homogeneous Neumann boundary conditions
    if solver.D1 isa DerivativeOperator ||
        solver.D1 isa UniformCoupledOperator
        D1_b = BandedMatrix(solver.D1)
        invImDKDn = inv(I + 1 / 6 * inv(M) * D1_b' * PdM * K * D1_b)
    elseif solver.D1 isa UpwindOperators
        D1plus_b = BandedMatrix(solver.D1.plus)
        invImDKDn = inv(I + 1 / 6 * inv(M) * D1plus_b' * PdM * K * D1plus_b)
    else
        @error "unknown type of first-derivative operator: $(typeof(solver.D1))"
    end
    return (invImD2Kd = invImD2Kd, invImDKDn = invImDKDn, D = D)
end

# Discretization that conserves the mass (for eta and v) and the energy for periodic boundary conditions, see
# - Joshua Lampert and Hendrik Ranocha (2024)
#   Structure-Preserving Numerical Methods for Two Nonlinear Systems of Dispersive Wave Equations
#   [DOI: 10.48550/arXiv.2402.16669](https://doi.org/10.48550/arXiv.2402.16669)
function rhs!(du_ode, u_ode, t, mesh, equations::BBMBBMVariableEquations1D,
              initial_condition, ::BoundaryConditionPeriodic, source_terms,
              solver, cache)
    @unpack invImDKD, invImD2K, D = cache

    q = wrap_array(u_ode, mesh, equations, solver)
    dq = wrap_array(du_ode, mesh, equations, solver)

    eta = view(q, 1, :)
    v = view(q, 2, :)
    deta = view(dq, 1, :)
    dv = view(dq, 2, :)
    dD = view(dq, 3, :)
    fill!(dD, zero(eltype(dD)))

    if solver.D1 isa PeriodicDerivativeOperator ||
       solver.D1 isa UniformPeriodicCoupledOperator
        @timeit timer() "deta hyperbolic" deta[:]=-solver.D1 * (D .* v + eta .* v)
        @timeit timer() "dv hyperbolic" dv[:]=-solver.D1 *
                                              (equations.gravity * eta + 0.5 * v .^ 2)
    elseif solver.D1 isa PeriodicUpwindOperators
        @timeit timer() "deta hyperbolic" deta[:]=-solver.D1.minus * (D .* v + eta .* v)
        @timeit timer() "dv hyperbolic" dv[:]=-solver.D1.plus *
                                              (equations.gravity * eta + 0.5 * v .^ 2)
    else
        @error "unknown type of first-derivative operator: $(typeof(solver.D1))"
    end

    @timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations, solver)

    @timeit timer() "deta elliptic" deta[:]=invImDKD * deta
    @timeit timer() "dv elliptic" dv[:]=invImD2K * dv

    return nothing
end

function rhs!(du_ode, u_ode, t, mesh, equations::BBMBBMVariableEquations1D,
              initial_condition, ::BoundaryConditionReflecting, source_terms,
              solver, cache)
    @unpack invImDKDn, invImD2Kd, D = cache

    q = wrap_array(u_ode, mesh, equations, solver)
    dq = wrap_array(du_ode, mesh, equations, solver)

    eta = view(q, 1, :)
    v = view(q, 2, :)
    deta = view(dq, 1, :)
    dv = view(dq, 2, :)
    dD = view(dq, 3, :)
    fill!(dD, zero(eltype(dD)))

    # energy and mass conservative semidiscretization
    if solver.D1 isa DerivativeOperator ||
       solver.D1 isa UniformCoupledOperator
        @timeit timer() "deta hyperbolic" deta[:]=-solver.D1 * (D .* v + eta .* v)
        @timeit timer() "dv hyperbolic" dv[:]=-solver.D1 *
                                            (equations.gravity * eta + 0.5 * v .^ 2)
    elseif solver.D1 isa UpwindOperators
        @timeit timer() "deta hyperbolic" deta[:]=-solver.D1.minus * (D .* v + eta .* v)
        @timeit timer() "dv hyperbolic" dv[:]=-solver.D1.plus *
                                            (equations.gravity * eta + 0.5 * v .^ 2)
    else
        @error "unknown type of first-derivative operator: $(typeof(solver.D1))"
    end

    @timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations, solver)

    @timeit timer() "deta elliptic" deta[:]=invImDKDn * deta
    @timeit timer() "dv elliptic" begin
    dv[2:(end - 1)] = invImD2Kd * dv[2:(end - 1)]
    dv[1] = dv[end] = zero(eltype(dv))
end

return nothing
end

@inline function prim2cons(q, equations::BBMBBMVariableEquations1D)
    eta, v, D = q

    b = bathymetry(q, equations)
    h = eta - b
    hv = h * v
    return SVector(h, hv, b)
end

@inline function cons2prim(u, equations::BBMBBMVariableEquations1D)
    h, hv, b = u

    eta = h + b
    v = hv / h
    D = equations.eta0 - b
    return SVector(eta, v, D)
end

@inline function waterheight_total(q, equations::BBMBBMVariableEquations1D)
    return q[1]
end

@inline function velocity(q, equations::BBMBBMVariableEquations1D)
    return q[2]
end

@inline function bathymetry(q, equations::BBMBBMVariableEquations1D)
    D = q[3]
    return equations.eta0 - D
end

@inline function waterheight(q, equations::BBMBBMVariableEquations1D)
    return waterheight_total(q, equations) - bathymetry(q, equations)
end

@inline function energy_total(q, equations::BBMBBMVariableEquations1D)
    eta, v, D = q
    e = 0.5 * (equations.gravity * eta^2 + (D + eta - equations.eta0) * v^2)
    return e
end

@inline entropy(q, equations::BBMBBMVariableEquations1D) = energy_total(q, equations)

# Calculate the error for the "lake-at-rest" test case where eta should
# be a constant value over time
@inline function lake_at_rest_error(q, equations::BBMBBMVariableEquations1D)
    eta, _, _ = q
    return abs(equations.eta0 - eta)
end
