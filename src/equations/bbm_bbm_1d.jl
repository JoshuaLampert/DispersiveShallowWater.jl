@doc raw"""
    BBMBBMEquations1D(; bathymetry_type = bathymetry_variable,
                      gravity, eta0 = 0.0)

BBM-BBM (Benjamin–Bona–Mahony) system in one spatial dimension. The equations for flat bathymetry are given by
```math
\begin{aligned}
  \eta_t + ((\eta + D)v)_x - \frac{1}{6}D^2\eta_{xxt} &= 0,\\
  v_t + g\eta_x + \left(\frac{1}{2}v^2\right)_x - \frac{1}{6}D^2v_{xxt} &= 0.
\end{aligned}
```
The unknown quantities of the BBM-BBM equations are the total water height ``\eta`` and the velocity ``v``.
The gravitational acceleration `gravity` is denoted by ``g`` and the constant bottom topography (bathymetry) ``b = \eta_0 - D``.
The water height above the bathymetry is therefore given by ``h = \eta - \eta_0 + D``.
The BBM-BBM equations are only implemented for ``\eta_0 = 0``.

Two types of `bathymetry_type` are supported:
- [`bathymetry_flat`](@ref): flat bathymetry (typically ``b = 0`` everywhere)
- [`bathymetry_variable`](@ref): general variable bathymetry

For the general case of variable vathymetry the BBM-BBM equations are
```math
\begin{aligned}
  \eta_t + ((\eta + D)v)_x - \frac{1}{6}(D^2\eta_{xt})_x &= 0,\\
  v_t + g\eta_x + \left(\frac{1}{2}v^2\right)_x - \frac{1}{6}(D^2v_t)_{xx} &= 0.
\end{aligned}
```

One reference for the BBM-BBM system can be found in Bona et al. (1998).
The semidiscretization implemented here was developed for flat bathymetry in
Ranocha et al. (2020) and generalized for a variable bathymetry in
Lampert and Ranocha (2024). It conserves
- the total water mass (integral of ``h``) as a linear invariant
- the total velocity (integral of ``v``) as a linear invariant for flat bathymetry
- the total energy

for periodic boundary conditions (see Lampert, Ranocha). For reflecting boundary conditions,
the semidiscretization conserves
- the total water (integral of ``h``) as a linear invariant
- the total energy.

Additionally, it is well-balanced for the lake-at-rest stationary solution, see Lampert and Ranocha (2024).

- Jerry L. Bona, Min Chen (1998)
  A Boussinesq system for two-way propagation of nonlinear dispersive waves
  [DOI: 10.1016/S0167-2789(97)00249-2](https://doi.org/10.1016/S0167-2789(97)00249-2)
- Hendrik Ranocha, Dimitrios Mitsotakis, David I. Ketcheson (2020)
  A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations
  [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)
- Joshua Lampert, Hendrik Ranocha (2024)
  Structure-Preserving Numerical Methods for Two Nonlinear Systems of Dispersive Wave Equations
  [DOI: 10.48550/arXiv.2402.16669](https://doi.org/10.48550/arXiv.2402.16669)
"""
struct BBMBBMEquations1D{Bathymetry <: AbstractBathymetry, RealT <: Real} <:
       AbstractBBMBBMEquations{1, 3}
    bathymetry_type::Bathymetry # type of bathymetry
    gravity::RealT # gravitational acceleration
    eta0::RealT    # constant still-water surface
end

function BBMBBMEquations1D(; bathymetry_type = bathymetry_variable,
                           gravity, eta0 = 0.0)
    eta0 == 0.0 || @warn "The still-water surface needs to be 0 for the BBM-BBM equations"
    BBMBBMEquations1D(bathymetry_type, gravity, eta0)
end

"""
    initial_condition_convergence_test(x, t, equations::BBMBBMEquations1D, mesh)

A traveling-wave solution used for convergence tests in a periodic domain.
The bathymetry is constant.

For details see Example 5 in Section 3 from (here adapted for dimensional equations):
- Min Chen (1997)
  Exact Traveling-Wave Solutions to Bidirectional Wave Equations
  [DOI: 10.1023/A:1026667903256](https://doi.org/10.1023/A:1026667903256)
"""
function initial_condition_convergence_test(x, t, equations::BBMBBMEquations1D, mesh)
    g = gravity(equations)
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
    initial_condition_manufactured(x, t, equations::BBMBBMEquations1D, mesh)

A smooth manufactured solution in combination with [`source_terms_manufactured`](@ref).
"""
function initial_condition_manufactured(x, t,
                                        equations::BBMBBMEquations1D,
                                        mesh)
    eta = exp(t) * cospi(2 * (x - 2 * t))
    v = exp(t / 2) * sinpi(2 * (x - t / 2))
    if equations.bathymetry_type isa BathymetryFlat
        D = 4.0
    else # equations.bathymetry_type isa BathymetryVariable
        D = 5 + 2 * cospi(2 * x)
    end
    return SVector(eta, v, D)
end

"""
    source_terms_manufactured(q, x, t, equations::BBMBBMEquations1D, mesh)

A smooth manufactured solution in combination with [`initial_condition_manufactured`](@ref).
"""
function source_terms_manufactured(q, x, t, equations::BBMBBMEquations1D)
    g = gravity(equations)
    a1 = cospi(2 * x)
    a2 = sinpi(2 * x)
    a3 = cospi(t - 2 * x)
    a4 = sinpi(t - 2 * x)
    a5 = sinpi(2 * t - 4 * x)
    a6 = sinpi(4 * t - 2 * x)
    a7 = cospi(4 * t - 2 * x)
    a9 = sinpi(x)
    a10 = exp(t)
    a14 = exp(t / 2)
    a15 = exp(3 * t / 2)
    dq1 = -2 * pi^2 * (4 * pi * a6 - a7) * (2 * a1 + 5)^2 * a10 / 3 +
          8 * pi^2 * (a6 + 4 * pi * a7) * (2 * a1 + 5) * a10 * a2 / 3 +
          2 * pi * (2 * a1 + 5) * a14 * a3 - 2 * pi * a15 * a4 * a6 +
          2 * pi * a15 * a3 * a7 + 4 * pi * a14 * a2 * a4 -
          4 * pi * a10 * a6 + a10 * a7
    dq2 = 2 * pi * g * a10 * a6 -
          pi^2 *
          (8 * (2 * pi * a4 - a3) * (2 * a1 + 5) * a2 +
           (a4 + 2 * pi * a3) * (2 * a1 + 5)^2 +
           4 * (a4 + 2 * pi * a3) * (16 * a9^4 - 26 * a9^2 + 7)) * a14 /
          3 - a14 * a4 / 2 - pi * a14 * a3 - pi * a10 * a5

    return SVector(dq1, dq2, zero(dq1))
end

function source_terms_manufactured(q, x, t, equations::BBMBBMEquations1D{BathymetryFlat})
    g = gravity(equations)
    D = still_waterdepth(q, equations)
    a3 = cospi(t - 2 * x)
    a4 = sinpi(t - 2 * x)
    a5 = sinpi(2 * t - 4 * x)
    a6 = sinpi(4 * t - 2 * x)
    a7 = cospi(4 * t - 2 * x)
    a10 = exp(t)
    a14 = exp(t / 2)
    a15 = exp(3 * t / 2)
    dq1 = -2 * pi^2 * D^2 * (4 * pi * a6 - a7) * a10 / 3 + 2 * pi * D * a14 * a3 -
          2 * pi * a15 * a4 * a6 + 2 * pi * a15 * a3 * a7 -
          4 * pi * a10 * a6 + a10 * a7
    dq2 = -pi^2 * D^2 * (a4 + 2 * pi * a3) * a14 / 3 + 2 * pi * g * a10 * a6 -
          a14 * a4 / 2 - pi * a14 * a3 - pi * a10 * a5

    return SVector(dq1, dq2, zero(dq1))
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
    if equations.bathymetry_type isa BathymetryFlat
        D = 4.0
    else # equations.bathymetry_type isa BathymetryVariable
        D = 5 + 2 * cospi(2 * x)
    end
    return SVector(eta, v, D)
end

"""
    source_terms_manufactured_reflecting(q, x, t, equations::BBMBBMEquations1D, mesh)

A smooth manufactured solution for reflecting boundary conditions in combination
with [`initial_condition_manufactured_reflecting`](@ref).
"""
function source_terms_manufactured_reflecting(q, x, t, equations::BBMBBMEquations1D)
    g = gravity(equations)
    a1 = cospi(2 * x)
    a2 = sinpi(2 * x)
    a8 = cospi(x)
    a9 = sinpi(x)
    a10 = exp(t)
    a11 = exp(2 * t)
    a12 = cospi(3 * x)
    a13 = sinpi(3 * x)
    dq1 = (pi * x * a11 * a1 + 4 * pi * x * a8 + 3 * pi * x * a12 +
           20 * pi^2 * (1 - a1)^2 * a10 * a8 / 3 + a11 * a2 / 2 + 2 * a10 * a8 +
           7 * pi^2 * a10 * a8 / 3 + 14 * pi^2 * a10 * a12 + 4 * a9 + a13) * a10
    dq2 = (-pi * g * a10 * a9 + pi * x^2 * a10 * a2 / 2 + x * a10 * a9^2 + x * a9 +
           pi * (400 * pi * x * a9^5 - 824 * pi * x * a9^3 + 385 * pi * x * a9 -
            160 * a9^4 * a8 + 336 * a9^2 * a8 - 98 * a8) / 6) * a10

    return SVector(dq1, dq2, zero(dq1))
end

function source_terms_manufactured_reflecting(q, x, t,
                                              equations::BBMBBMEquations1D{BathymetryFlat})
    g = gravity(equations)
    D = still_waterdepth(q, equations)
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

    return SVector(dq1, dq2, zero(dq1))
end

"""
    initial_condition_dingemans(x, t, equations::BBMBBMEquations1D, mesh)

The initial condition that uses the dispersion relation of the Euler equations
to approximate waves generated by a wave maker as it is done by experiments of
Dingemans. The topography is a trapezoidal.

!!! warning "Translation of water height"
    The initial condition for the water height is translated to be around 0, which is
    needed for the simulation because the `BBMBBMEquations1D` are only implemented
    for ``\\eta_0 = 0``.

References:
- Magnus Svärd, Henrik Kalisch (2023)
  A novel energy-bounded Boussinesq model and a well-balanced and stable numerical discretization
  [arXiv: 2302.09924](https://arxiv.org/abs/2302.09924)
- Maarten W. Dingemans (1994)
  Comparison of computations with Boussinesq-like models and laboratory measurements
  [link](https://repository.tudelft.nl/islandora/object/uuid:c2091d53-f455-48af-a84b-ac86680455e9/datastream/OBJ/download)
"""
function initial_condition_dingemans(x, t, equations::BBMBBMEquations1D, mesh)
    g = gravity(equations)
    h0 = 0.8
    A = 0.02
    # omega = 2*pi/(2.02*sqrt(2))
    k = 0.8406220896381442 # precomputed result of find_zero(k -> omega^2 - g * k * tanh(k * h0), 1.0) using Roots.jl
    if x < -30.5 * pi / k || x > -8.5 * pi / k
        h = 0.0
    else
        h = A * cos(k * x)
    end
    v = sqrt(g / k * tanh(k * h0)) * h / h0
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

function create_cache(mesh, equations::BBMBBMEquations1D{BathymetryFlat},
                      solver, initial_condition,
                      ::BoundaryConditionPeriodic,
                      RealT, uEltype)
    # D is constant for flat bathymetry
    x = grid(solver)
    D = still_waterdepth(initial_condition(first(x), 0.0, equations, mesh), equations)
    invImD2 = lu(I - 1 / 6 * D^2 * sparse(solver.D2))

    # create temporary storage
    etav = zeros(RealT, nnodes(mesh))
    Dv = zero(etav)
    v2 = zero(etav)
    tmp2 = zero(etav)
    return (; invImD2, etav, Dv, v2, tmp2)
end

function create_cache(mesh, equations::BBMBBMEquations1D,
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
       solver.D1 isa UniformPeriodicCoupledOperator ||
       solver.D1 isa FourierDerivativeOperator
        sparse_D1 = sparse(solver.D1)
        invImDKD = lu(I - 1 / 6 * sparse_D1 * K * sparse_D1)
    elseif solver.D1 isa PeriodicUpwindOperators
        invImDKD = lu(I - 1 / 6 * sparse(solver.D1.minus) * K * sparse(solver.D1.plus))
    else
        throw(ArgumentError("unknown type of first-derivative operator: $(typeof(solver.D1))"))
    end
    invImD2K = lu(I - 1 / 6 * sparse(solver.D2) * K)

    # create temporary storage
    etav = zeros(RealT, nnodes(mesh))
    Dv = zero(etav)
    v2 = zero(etav)
    tmp2 = zero(etav)
    return (; invImDKD, invImD2K, etav, Dv, v2, tmp2)
end

function create_cache(mesh, equations::BBMBBMEquations1D{BathymetryFlat},
                      solver, initial_condition,
                      ::BoundaryConditionReflecting,
                      RealT, uEltype)
    # D is constant for flat bathymetry
    x = grid(solver)
    D = still_waterdepth(initial_condition(first(x), 0.0, equations, mesh), equations)
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
        throw(ArgumentError("unknown type of first-derivative operator: $(typeof(solver.D1))"))
    end

    # create temporary storage
    etav = zeros(RealT, nnodes(mesh))
    Dv = zero(etav)
    v2 = zero(etav)
    tmp2 = zero(etav)
    return (; invImD2d, invImD2n, etav, Dv, v2, tmp2)
end

function create_cache(mesh, equations::BBMBBMEquations1D,
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
    D2d = (BandedMatrix(solver.D2) * Pd)[2:(end - 1), :]
    # homogeneous Dirichlet boundary conditions
    invImD2d = lu(I - 1 / 6 * D2d * K_i)
    m = diag(M)
    m[1] = 0
    m[end] = 0
    PdM = Diagonal(m)

    # homogeneous Neumann boundary conditions
    if solver.D1 isa DerivativeOperator ||
       solver.D1 isa UniformCoupledOperator
        D1_b = BandedMatrix(solver.D1)
        invImD2n = lu(I + 1 / 6 * inv(M) * D1_b' * PdM * K * D1_b)
    elseif solver.D1 isa UpwindOperators
        D1plus_b = BandedMatrix(solver.D1.plus)
        invImD2n = lu(I + 1 / 6 * inv(M) * D1plus_b' * PdM * K * D1plus_b)
    else
        throw(ArgumentError("unknown type of first-derivative operator: $(typeof(solver.D1))"))
    end

    # create temporary storage
    etav = zeros(RealT, nnodes(mesh))
    Dv = zero(etav)
    v2 = zero(etav)
    tmp2 = zero(etav)
    return (; invImD2d, invImD2n, etav, Dv, v2, tmp2)
end

# Discretization that conserves
# - the total water (integral of ``h``) as a linear invariant
# - the total momentum (integral of ``v``) as a linear invariant for flat bathymetry
# - the total energy
# for periodic boundary conditions, see
# - Joshua Lampert and Hendrik Ranocha (2024)
#   Structure-Preserving Numerical Methods for Two Nonlinear Systems of Dispersive Wave Equations
#   [DOI: 10.48550/arXiv.2402.16669](https://doi.org/10.48550/arXiv.2402.16669)
function rhs!(dq, q, t, mesh, equations::BBMBBMEquations1D, initial_condition,
              ::BoundaryConditionPeriodic, source_terms, solver, cache)
    (; etav, Dv, v2, tmp1, tmp2) = cache
    if equations.bathymetry_type isa BathymetryFlat
        (; invImD2) = cache
    else # equations.bathymetry_type isa BathymetryVariable
        (; invImDKD, invImD2K) = cache
    end

    g = gravity(equations)
    eta, v, D = q.x
    deta, dv, dD = dq.x
    fill!(dD, zero(eltype(dD)))

    @trixi_timeit timer() "deta hyperbolic" begin
        @.. etav = eta * v
        @.. Dv = D * v
        @.. tmp1 = -Dv - etav
    end
    @trixi_timeit timer() "dv hyperbolic" begin
        @.. v2 = v .^ 2
        @.. tmp2 = -g * eta - 0.5 * v2
    end
    # energy and mass conservative semidiscretization
    if solver.D1 isa PeriodicDerivativeOperator ||
       solver.D1 isa UniformPeriodicCoupledOperator ||
       solver.D1 isa FourierDerivativeOperator
        @trixi_timeit timer() "deta hyperbolic" begin
            mul!(deta, solver.D1, tmp1)
        end
        @trixi_timeit timer() "dv hyperbolic" begin
            mul!(dv, solver.D1, tmp2)
        end
    elseif solver.D1 isa PeriodicUpwindOperators
        @trixi_timeit timer() "deta hyperbolic" begin
            mul!(deta, solver.D1.minus, tmp1)
        end
        @trixi_timeit timer() "dv hyperbolic" begin
            mul!(dv, solver.D1.plus, tmp2)
        end
    else
        throw(ArgumentError("unknown type of first-derivative operator: $(typeof(solver.D1))"))
    end

    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)

    @trixi_timeit timer() "solving elliptic system deta" begin
        if equations.bathymetry_type isa BathymetryFlat
            solve_system_matrix!(deta, invImD2, equations)
        else # equations.bathymetry_type isa BathymetryVariable
            solve_system_matrix!(deta, invImDKD, equations)
        end
    end
    @trixi_timeit timer() "solving elliptic system dv" begin
        if equations.bathymetry_type isa BathymetryFlat
            solve_system_matrix!(dv, invImD2, equations)
        else # equations.bathymetry_type isa BathymetryVariable
            solve_system_matrix!(dv, invImD2K, equations)
        end
    end
    return nothing
end

# Discretization that conserves
# - the total water (integral of ``h``) as a linear invariant
# - the total energy
# for reflecting boundary conditions, see
# - Joshua Lampert and Hendrik Ranocha (2024)
#   Structure-Preserving Numerical Methods for Two Nonlinear Systems of Dispersive Wave Equations
#   [DOI: 10.48550/arXiv.2402.16669](https://doi.org/10.48550/arXiv.2402.16669)
function rhs!(dq, q, t, mesh, equations::BBMBBMEquations1D, initial_condition,
              ::BoundaryConditionReflecting, source_terms, solver, cache)
    (; etav, Dv, v2, tmp1, tmp2) = cache
    (; invImD2d, invImD2n) = cache

    g = gravity(equations)
    eta, v, D = q.x
    deta, dv, dD = dq.x
    fill!(dD, zero(eltype(dD)))

    @trixi_timeit timer() "deta hyperbolic" begin
        @.. etav = eta * v
        @.. Dv = D * v
        @.. tmp1 = -Dv - etav
    end
    @trixi_timeit timer() "dv hyperbolic" begin
        @.. v2 = v .^ 2
        @.. tmp2 = -g * eta - 0.5 * v2
    end
    # energy and mass conservative semidiscretization
    if solver.D1 isa DerivativeOperator ||
       solver.D1 isa UniformCoupledOperator
        @trixi_timeit timer() "deta hyperbolic" begin
            mul!(deta, solver.D1, tmp1)
        end
        @trixi_timeit timer() "dv hyperbolic" begin
            mul!(dv, solver.D1, tmp2)
        end
    elseif solver.D1 isa UpwindOperators
        @trixi_timeit timer() "deta hyperbolic" begin
            mul!(deta, solver.D1.minus, tmp1)
        end
        @trixi_timeit timer() "dv hyperbolic" begin
            mul!(dv, solver.D1.plus, tmp2)
        end
    else
        throw(ArgumentError("unknown type of first-derivative operator: $(typeof(solver.D1))"))
    end

    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)

    @trixi_timeit timer() "solving elliptic system deta" begin
        solve_system_matrix!(deta, invImD2n, equations)
    end
    @trixi_timeit timer() "solving elliptic system dv" begin
        solve_system_matrix!((@view dv[2:(end - 1)]), invImD2d, equations)
        dv[1] = dv[end] = zero(eltype(dv))
    end

    return nothing
end
