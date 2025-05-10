@doc raw"""
    SvaerdKalischEquations1D(; bathymetry_type = bathymetry_variable, gravity,
                             eta0 = 0.0, alpha = 0.0, beta = 1/3, gamma = 0.0)

Dispersive system by Svärd and Kalisch (2023) in one spatial dimension. The equations for variable bathymetry
are given in conservative variables by
```math
\begin{aligned}
  h_t + (hv)_x &= (\hat\alpha(\hat\alpha(h + b)_x)_x)_x,\\
  (hv)_t + (hv^2)_x + gh(h + b)_x &= (\hat\alpha v(\hat\alpha(h + b)_x)_x)_x + (\hat\beta v_x)_{xt} + \frac{1}{2}(\hat\gamma v_x)_{xx} + \frac{1}{2}(\hat\gamma v_{xx})_x,
\end{aligned}
```
where ``\hat\alpha^2 = \alpha\sqrt{gD}D^2``, ``\hat\beta = \beta D^3``, ``\hat\gamma = \gamma\sqrt{gD}D^3``.
The coefficients ``\alpha``, ``\beta`` and ``\gamma`` are provided in dimensionless form and ``D = \eta_0 - b`` is the
still-water depth and `eta0` is the still-water surface (lake-at-rest).
The equations can be rewritten in primitive variables as
```math
\begin{aligned}
  \eta_t + ((\eta + D)v)_x &= (\hat\alpha(\hat\alpha\eta_x)_x)_x,\\
  v_t(\eta + D) - v((\eta + D)v)_x + ((\eta + D)v^2)_x + g(\eta + D)\eta_x &= (\hat\alpha v(\hat\alpha\eta_x)_x)_x - v(\hat\alpha(\hat\alpha\eta_x)_x)_x + (\hat\beta v_x)_{xt} + \frac{1}{2}(\hat\gamma v_x)_{xx} + \frac{1}{2}(\hat\gamma v_{xx})_x.
\end{aligned}
```
The unknown quantities of the Svärd-Kalisch equations are the total water height ``\eta`` and the velocity ``v``.
The gravitational acceleration `gravity` is denoted by ``g`` and the bottom topography (bathymetry) ``b = \eta_0 - D``.
The water height above the bathymetry is therefore given by
``h = \eta - \eta_0 + D``.

Currently, the equations only support a general variable bathymetry, see [`bathymetry_variable`](@ref).

`SvärdKalischEquations1D` is an alias for `SvaerdKalischEquations1D`.

The equations by Svärd and Kalisch are presented and analyzed in Svärd and Kalisch (2025).
The semidiscretization implemented here conserves
- the total water mass (integral of ``h``) as a linear invariant
- the total momentum (integral of ``h v``) as a nonlinear invariant for flat bathymetry
- the total modified energy

for periodic boundary conditions (see Lampert, Ranocha).
Additionally, it is well-balanced for the lake-at-rest stationary solution, see Lampert and Ranocha (2024).

- Magnus Svärd, Henrik Kalisch (2025)
  A novel energy-bounded Boussinesq model and a well-balanced and stable numerical discretization
  [arXiv: 2302.09924](https://arxiv.org/abs/2302.09924),
  [DOI: 10.1016/j.jcp.2024.113516](https://doi.org/10.1016/j.jcp.2024.113516)
- Joshua Lampert, Hendrik Ranocha (2024)
  Structure-Preserving Numerical Methods for Two Nonlinear Systems of Dispersive Wave Equations
  [DOI: 10.48550/arXiv.2402.16669](https://doi.org/10.48550/arXiv.2402.16669)
"""
struct SvaerdKalischEquations1D{Bathymetry <: AbstractBathymetry, RealT <: Real} <:
       AbstractSvaerdKalischEquations{1, 3}
    bathymetry_type::Bathymetry # type of bathymetry
    gravity::RealT # gravitational acceleration
    eta0::RealT    # constant still-water surface
    alpha::RealT   # coefficient
    beta::RealT    # coefficient
    gamma::RealT   # coefficient
end

"""
    SvärdKalischEquations1D

Same as [`SvaerdKalischEquations1D`](@ref).
"""
const SvärdKalischEquations1D = SvaerdKalischEquations1D

function SvaerdKalischEquations1D(; bathymetry_type = bathymetry_variable, gravity,
                                  eta0 = 0.0, alpha = 0.0, beta = 1 / 3, gamma = 0.0)
    SvaerdKalischEquations1D(bathymetry_type, gravity, eta0, alpha, beta, gamma)
end

"""
    initial_condition_manufactured(x, t, equations::SvaerdKalischEquations1D, mesh)

A smooth manufactured solution in combination with [`source_terms_manufactured`](@ref).
"""
function initial_condition_manufactured(x, t,
                                        equations::SvaerdKalischEquations1D,
                                        mesh)
    eta = exp(t) * cospi(2 * (x - 2 * t))
    v = exp(t / 2) * sinpi(2 * (x - t / 2))
    b = -5 - 2 * cospi(2 * x)
    D = equations.eta0 - b
    return SVector(eta, v, D)
end

"""
    source_terms_manufactured(q, x, t, equations::SvaerdKalischEquations1D, mesh)

A smooth manufactured solution in combination with [`initial_condition_manufactured`](@ref).
"""
function source_terms_manufactured(q, x, t, equations::SvaerdKalischEquations1D)
    g = gravity(equations)
    eta0 = still_water_surface(q, equations)
    alpha = equations.alpha
    beta = equations.beta
    gamma = equations.gamma
    a1 = sinpi(2 * x)
    a2 = cospi(2 * x)
    a3 = sinpi(-t + 2 * x)
    a4 = cospi(-t + 2 * x)
    a5 = sinpi(t - 2 * x)
    a6 = cospi(t - 2 * x)
    a7 = sinpi(-4 * t + 2 * x)
    a8 = exp(t / 2)
    a9 = exp(t)
    a10 = a9 * cospi(-4 * t + 2 * x)
    a11 = eta0 + 2.0 * a2 + 5.0
    a12 = sqrt(g * a11)
    a13 = 0.2 * eta0 + 0.4 * a2 + 1
    a14 = alpha * a12 * a13^2
    a15 = sqrt(a14)
    a16 = -1.0 * pi * a14 * a1 / a11 - 0.8 * pi * alpha * a12 * a13 * a1
    a17 = -20.0 * pi^2 * a15 * a10 - 10.0 * pi * a15 * a16 * a9 * a7 / (a14)
    a18 = -2 * pi * a9 * a7 - 4.0 * pi * a1
    a19 = a10 + 2.0 * a2 + 5.0
    a20 = a18 * a8 * a3 + 2 * pi * a19 * a8 * a4
    a21 = a15 * (40.0 * pi^3 * a15 * a9 * a7 - 40.0 * pi^2 * a15 * a16 * a10 / (a14) -
           20.0 * pi^2 * a15 * a16 * a9 * a1 * a7 / (a14 * a11) -
           16.0 * pi^2 * a15 * a16 * a9 * a1 * a7 / (alpha * a12 * a13^3) -
           10.0 * pi * a15 *
           (-2.0 * pi^2 * a14 * a2 / a11 - 1.6 * pi^2 * alpha * a12 * a13 * a2 +
            3.2 * pi^2 * alpha * a12 * a13 * a1^2 / a11 + 0.56 * pi^2 * alpha * a12 * a1^2) *
           a9 * a7 / (a14) -
           10.0 * pi * a15 * a16^2 * a9 * a7 / (alpha^2 * g * a13^4 * a11))

    dq1 = -5.0 * a21 + a20 + 4 * pi * a9 * a7 + a10 - 5.0 * a15 * a17 * a16 / (a14)

    dq2 = -25.0 * beta * (-2 * pi^2 * a8 * a3 + 4 * pi^3 * a8 * a4) * a13^2 * a11 +
          100.0 * pi * beta * (2 * pi^2 * a8 * a3 + pi * a8 * a4) * a13^2 * a1 +
          40.0 * pi * beta * (2 * pi^2 * a8 * a3 + pi * a8 * a4) * a13 * a11 * a1 -
          2 * pi * g * a19 * a9 * a7 + 100.0 * pi^3 * gamma * a12 * a13^2 * a11 * a8 * a4 -
          300.0 * pi^3 * gamma * a12 * a13^2 * a8 * a1 * a3 -
          80.0 * pi^3 * gamma * a12 * a13 * a11 * a8 * a1 * a3 -
          pi^3 * gamma * a12 *
          (-50.0 * (3.2 * a13 * a2 - 1.28 * a1^2) * a11 * a6 -
           50.0 * (4.0 * a2 / a11 + 0.16 * a1^2 / a13^2) * a13^2 * a11 * a6 -
           200.0 * a13^2 * a11 * a6 - 1200.0 * a13^2 * a1 * a5 - 400.0 * a13^2 * a2 * a6 +
           800.0 * a13^2 * a1^2 * a6 / a11 - 320.0 * a13 * a11 * a1 * a5 +
           960.0 * a13 * a1^2 * a6) * a8 / 2 - 10.0 * pi * a15 * a17 * a8 * a4 -
          2.5 * a21 * a8 * a3 + (5.0 * a21 + 5.0 * a15 * a17 * a16 / (a14)) * a8 * a3 / 2 +
          (a8 * a3 / 2 - pi * a8 * a4) * a19 + a18 * a9 * a3^2 / 2 - (a20) * a8 * a3 / 2 +
          3 * pi * a19 * a9 * a3 * a4 - 2.5 * a15 * a17 * a16 * a8 * a3 / (a14)

    return SVector(dq1, dq2, zero(dq1))
end

"""
    initial_condition_manufactured_reflecting(x, t, equations::SvaerdKalischEquations1D, mesh)

A smooth manufactured solution for reflecting boundary conditions in combination
with [`source_terms_manufactured_reflecting`](@ref).
"""
function initial_condition_manufactured_reflecting(x, t,
                                                   equations::SvaerdKalischEquations1D,
                                                   mesh)
    eta = exp(2 * t) * cospi(x)
    v = exp(t) * x * sinpi(x)
    b = -5 - 2 * cospi(2 * x)
    D = equations.eta0 - b
    return SVector(eta, v, D)
end

"""
    source_terms_manufactured_reflecting(q, x, t, equations::SvaerdKalischEquations1D, mesh)

A smooth manufactured solution for reflecting boundary conditions in combination
with [`initial_condition_manufactured_reflecting`](@ref).
"""
function source_terms_manufactured_reflecting(q, x, t, equations::SvaerdKalischEquations1D)
    g = gravity(equations)
    eta0 = still_water_surface(q, equations)
    beta = equations.beta
    a1 = sinpi(2 * x)
    a2 = cospi(2 * x)
    a9 = exp(t)
    a11 = eta0 + 2.0 * a2 + 5.0
    a13 = 0.2 * eta0 + 0.4 * a2 + 1
    a22 = sinpi(x)
    a23 = cospi(x)
    a24 = exp(2 * t)
    a25 = a24 * a23
    a26 = a25 + 2.0 * a2 + 5.0
    a27 = a9 * a22
    a28 = pi * x * a9 * a23 + a27
    a29 = -pi * a24 * a22 - 4.0 * pi * a1
    a30 = a26 * a24
    a31 = a26 * a27
    dq1 = x * a29 * a27 + pi * x * a26 * a9 * a23 + a31 + 2 * a25
    dq2 = 100.0 * pi * beta * a28 * a13^2 * a1 + 40.0 * pi * beta * a28 * a13 * a11 * a1 -
          25.0 * beta * (-pi^2 * x * a27 + 2 * pi * a9 * a23) * a13^2 * a11 -
          pi * g * a30 * a22 + x^2 * a29 * a24 * a22^2 / 2 + pi * x^2 * a30 * a22 * a23 +
          x * a28 * a31 / 2 + x * a30 * a22^2 + x * a31 -
          x * (x * a29 * a27 + pi * x * a26 * a9 * a23 + a31) * a27 / 2

    return SVector(dq1, dq2, zero(dq1))
end

# For periodic boundary conditions
function assemble_system_matrix!(cache, h,
                                 ::SvaerdKalischEquations1D,
                                 ::BoundaryConditionPeriodic)
    (; D1, M_h, minus_MD1betaD1) = cache

    @.. M_h = h
    scale_by_mass_matrix!(M_h, D1)

    # Floating point errors accumulate a bit and the system matrix
    # is not necessarily perfectly symmetric but only up to
    # round-off errors. We wrap it here to avoid issues with the
    # factorization.
    return Symmetric(Diagonal(M_h) + minus_MD1betaD1)
end

# For reflecting boundary conditions
function assemble_system_matrix!(cache, h,
                                 ::SvaerdKalischEquations1D,
                                 ::BoundaryConditionReflecting)
    (; D1, M_h, minus_MD1betaD1) = cache

    @.. M_h = h
    scale_by_mass_matrix!(M_h, D1)

    # Floating point errors accumulate a bit and the system matrix
    # is not necessarily perfectly symmetric but only up to
    # round-off errors. We wrap it here to avoid issues with the
    # factorization.
    return Symmetric(Diagonal((@view M_h[(begin + 1):(end - 1)])) + minus_MD1betaD1)
end

function create_cache(mesh, equations::SvaerdKalischEquations1D,
                      solver, initial_condition,
                      boundary_conditions::BoundaryConditionPeriodic,
                      RealT, uEltype)
    D1 = solver.D1
    #  Assume D is independent of time and compute D evaluated at mesh points once.
    D = Array{RealT}(undef, nnodes(mesh))
    x = grid(solver)
    for i in eachnode(solver)
        D[i] = still_waterdepth(initial_condition(x[i], 0.0, equations, mesh), equations)
    end
    g = gravity(equations)
    h = ones(RealT, nnodes(mesh))
    hv = zero(h)
    b = zero(h)
    eta_x = zero(h)
    v_x = zero(h)
    alpha_eta_x_x = zero(h)
    y_x = zero(h)
    v_y_x = zero(h)
    vy = zero(h)
    vy_x = zero(h)
    y_v_x = zero(h)
    h_v_x = zero(h)
    hv2_x = zero(h)
    v_xx = zero(h)
    gamma_v_xx_x = zero(h)
    gamma_v_x_xx = zero(h)
    alpha_hat = sqrt.(equations.alpha * sqrt.(g * D) .* D .^ 2)
    beta_hat = equations.beta * D .^ 3
    gamma_hat = equations.gamma * sqrt.(g * D) .* D .^ 3
    tmp2 = zero(h)
    M_h = copy(h)
    scale_by_mass_matrix!(M_h, D1)
    M_beta = copy(beta_hat)
    scale_by_mass_matrix!(M_beta, D1)
    if D1 isa PeriodicDerivativeOperator ||
       D1 isa UniformPeriodicCoupledOperator ||
       D1 isa FourierDerivativeOperator
        D1_central = D1
        D1mat = sparse(D1_central)
        minus_MD1betaD1 = D1mat' * (Diagonal(M_beta) * D1mat)
        system_matrix = let cache = (; D1, M_h, minus_MD1betaD1)
            assemble_system_matrix!(cache, h,
                                    equations, boundary_conditions)
        end
    elseif D1 isa PeriodicUpwindOperators
        D1_central = D1.central
        D1mat_minus = sparse(D1.minus)
        minus_MD1betaD1 = D1mat_minus' * (Diagonal(M_beta) * D1mat_minus)
        system_matrix = let cache = (; D1, M_h, minus_MD1betaD1)
            assemble_system_matrix!(cache, h,
                                    equations, boundary_conditions)
        end
    else
        throw(ArgumentError("unknown type of first-derivative operator: $(typeof(D1))"))
    end
    factorization = cholesky(system_matrix)
    cache = (; factorization, minus_MD1betaD1, D, h, hv, b, eta_x, v_x,
             alpha_eta_x_x, y_x, v_y_x, vy, vy_x, y_v_x, h_v_x, hv2_x, v_xx,
             gamma_v_xx_x, gamma_v_x_xx,
             alpha_hat, beta_hat, gamma_hat, tmp2, D1_central, D1, M_h)
    if D1 isa PeriodicUpwindOperators
        eta_x_upwind = zero(h)
        v_x_upwind = zero(h)
        cache = (; cache..., eta_x_upwind, v_x_upwind)
    end
    return cache
end

# Reflecting boundary conditions assume alpha = gamma = 0
function create_cache(mesh, equations::SvaerdKalischEquations1D,
                      solver, initial_condition,
                      boundary_conditions::BoundaryConditionReflecting,
                      RealT, uEltype)
    if !iszero(equations.alpha) || !iszero(equations.gamma)
        throw(ArgumentError("Reflecting boundary conditions for Svärd-Kalisch equations only implemented for alpha = gamma = 0"))
    end
    D1 = solver.D1
    N = nnodes(mesh)
    #  Assume D is independent of time and compute D evaluated at mesh points once.
    D = Array{RealT}(undef, N)
    x = grid(solver)
    for i in eachnode(solver)
        D[i] = still_waterdepth(initial_condition(x[i], 0.0, equations, mesh), equations)
    end
    h = ones(RealT, N)
    hv = zero(h)
    b = zero(h)
    eta_x = zero(h)
    v_x = zero(h)
    h_v_x = zero(h)
    hv2_x = zero(h)
    beta_hat = equations.beta .* D .^ 3
    M_h = copy(h)
    scale_by_mass_matrix!(M_h, D1)
    M_beta = copy(beta_hat)
    scale_by_mass_matrix!(M_beta, D1)
    Pd = BandedMatrix((-1 => fill(one(real(mesh)), N - 2),), (N, N - 2))
    if D1 isa DerivativeOperator ||
       D1 isa UniformCoupledOperator
        D1_central = D1
        D1mat = sparse(D1_central)
        minus_MD1betaD1 = sparse(D1mat' * (Diagonal(M_beta) *
                                           D1mat * Pd))[(begin + 1):(end - 1), :]
        system_matrix = let cache = (; D1, M_h, minus_MD1betaD1)
            assemble_system_matrix!(cache, h, equations, boundary_conditions)
        end
    elseif D1 isa UpwindOperators
        D1_central = D1.central
        D1mat_minus = sparse(D1.minus)
        minus_MD1betaD1 = sparse(D1mat_minus' * (Diagonal(M_beta) *
                                                 D1mat_minus * Pd))[(begin + 1):(end - 1),
                                                                    :]
        system_matrix = let cache = (; D1, M_h, minus_MD1betaD1)
            assemble_system_matrix!(cache, h, equations, boundary_conditions)
        end
    else
        throw(ArgumentError("unknown type of first-derivative operator: $(typeof(D1))"))
    end
    factorization = cholesky(system_matrix)
    cache = (; factorization, minus_MD1betaD1, D, h, hv, b, eta_x, v_x,
             h_v_x, hv2_x, beta_hat, D1_central, D1, M_h)
    return cache
end

# Discretization that conserves
# - the total water (integral of ``h``) as a linear invariant
# - the total momentum (integral of ``h v``) as a nonlinear invariant for flat bathymetry
# - the total modified energy
# for periodic boundary conditions, see
# - Joshua Lampert and Hendrik Ranocha (2024)
#   Structure-Preserving Numerical Methods for Two Nonlinear Systems of Dispersive Wave Equations
#   [DOI: 10.48550/arXiv.2402.16669](https://doi.org/10.48550/arXiv.2402.16669)
# TODO: Simplify for the case of flat bathymetry and use higher-order operators
function rhs!(dq, q, t, mesh, equations::SvaerdKalischEquations1D,
              initial_condition, boundary_conditions::BoundaryConditionPeriodic,
              source_terms, solver, cache)
    (; D, h, hv, b, eta_x, v_x, alpha_eta_x_x, y_x, v_y_x, vy, vy_x,
    y_v_x, h_v_x, hv2_x, v_xx, gamma_v_xx_x, gamma_v_x_xx, alpha_hat, gamma_hat,
    tmp1, tmp2, D1_central, D1) = cache

    g = gravity(equations)
    eta, v = q.x
    deta, dv, dD = dq.x
    fill!(dD, zero(eltype(dD)))

    @trixi_timeit timer() "deta hyperbolic" begin
        @.. b = equations.eta0 - D
        @.. h = eta - b
        @.. hv = h * v

        if D1 isa PeriodicDerivativeOperator ||
           D1 isa UniformPeriodicCoupledOperator ||
           D1 isa FourierDerivativeOperator
            mul!(eta_x, D1_central, eta)
            mul!(v_x, D1_central, v)
            @.. tmp1 = alpha_hat * eta_x
            mul!(alpha_eta_x_x, D1_central, tmp1)
            @.. tmp1 = alpha_hat * alpha_eta_x_x
            mul!(y_x, D1_central, tmp1)
            @.. v_y_x = v * y_x
            @.. vy = v * tmp1
            mul!(vy_x, D1_central, vy)
            @.. y_v_x = tmp1 * v_x
            @.. tmp2 = tmp1 - hv
            mul!(deta, D1_central, tmp2)
        elseif D1 isa PeriodicUpwindOperators
            (; eta_x_upwind, v_x_upwind) = cache
            mul!(eta_x, D1_central, eta)
            mul!(v_x, D1_central, v)
            mul!(eta_x_upwind, D1.plus, eta)
            @.. tmp1 = alpha_hat * eta_x_upwind
            mul!(alpha_eta_x_x, D1.minus, tmp1)
            @.. tmp1 = alpha_hat * alpha_eta_x_x
            mul!(y_x, D1.minus, tmp1)
            @.. v_y_x = v * y_x
            @.. vy = v * tmp1
            mul!(vy_x, D1.minus, vy)
            mul!(v_x_upwind, D1.plus, v)
            @.. y_v_x = tmp1 * v_x_upwind
            # deta[:] = D1.minus * tmp1 - D1_central * hv
            mul!(deta, D1.minus, tmp1)
            mul!(deta, D1_central, hv, -1.0, 1.0)
        else
            throw(ArgumentError("unknown type of first-derivative operator: $(typeof(D1))"))
        end
    end

    # split form
    @trixi_timeit timer() "dv hyperbolic" begin
        mul!(h_v_x, D1_central, hv)
        @.. tmp1 = hv * v
        mul!(hv2_x, D1_central, tmp1)
        mul!(v_xx, solver.D2, v)
        @.. tmp1 = gamma_hat * v_xx
        mul!(gamma_v_xx_x, D1_central, tmp1)
        @.. tmp1 = gamma_hat * v_x
        mul!(gamma_v_x_xx, solver.D2, tmp1)
        @.. dv = -(0.5 * (hv2_x + hv * v_x - v * h_v_x) +
                   g * h * eta_x +
                   0.5 * (v_y_x - vy_x - y_v_x) -
                   0.5 * gamma_v_xx_x -
                   0.5 * gamma_v_x_xx)
    end

    # no split form
    #     dv[:] = -(D1_central * (hv .* v) - v .* (D1_central * hv)+
    #               g * h .* eta_x +
    #               vy_x - v_y_c -
    #               0.5 * gamma_v_xx_x -
    #               0.5 * gamma_v_x_xx)

    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)
    @trixi_timeit timer() "assemble system matrix" begin
        system_matrix = assemble_system_matrix!(cache, h, equations,
                                                boundary_conditions)
    end
    @trixi_timeit timer() "solving elliptic system" begin
        tmp1 .= dv
        solve_system_matrix!(dv, system_matrix,
                             tmp1, equations, D1, cache, boundary_conditions)
    end

    return nothing
end

function rhs!(dq, q, t, mesh, equations::SvaerdKalischEquations1D,
              initial_condition, boundary_conditions::BoundaryConditionReflecting,
              source_terms, solver, cache)
    # When constructing the `cache`, we assert that alpha and gamma are zero.
    # We use this explicitly in the code below.
    (; D, h, hv, b, eta_x, v_x, h_v_x, hv2_x, tmp1, D1_central, D1) = cache

    g = gravity(equations)
    eta, v = q.x
    deta, dv, dD = dq.x
    fill!(dD, zero(eltype(dD)))

    @trixi_timeit timer() "deta hyperbolic" begin
        @.. b = equations.eta0 - D
        @.. h = eta - b
        @.. hv = h * v

        mul!(deta, D1_central, -hv)
    end

    # split form
    @trixi_timeit timer() "dv hyperbolic" begin
        mul!(eta_x, D1_central, eta)
        mul!(v_x, D1_central, v)
        mul!(h_v_x, D1_central, hv)
        @.. tmp1 = hv * v
        mul!(hv2_x, D1_central, tmp1)
        @.. dv = -(0.5 * (hv2_x + hv * v_x - v * h_v_x) +
                   g * h * eta_x)
    end

    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)
    @trixi_timeit timer() "assemble system matrix" begin
        system_matrix = assemble_system_matrix!(cache, h, equations, boundary_conditions)
    end
    @trixi_timeit timer() "solving elliptic system" begin
        tmp1 .= dv
        solve_system_matrix!((@view dv[(begin + 1):(end - 1)]), system_matrix,
                             tmp1, equations, D1, cache, boundary_conditions)
        dv[1] = dv[end] = zero(eltype(dv))
    end

    return nothing
end

# The modified entropy/energy takes the whole `q` for every point in space
"""
    DispersiveShallowWater.energy_total_modified!(e, q_global, equations::SvaerdKalischEquations1D, cache)

Return the modified total energy `e` of the primitive variables `q_global` for the
[`SvaerdKalischEquations1D`](@ref). It contains an additional term containing a
derivative compared to the usual [`energy_total`](@ref) modeling
non-hydrostatic contributions. The `energy_total_modified`
is a conserved quantity (for periodic boundary conditions).

It is given by
```math
\\frac{1}{2} g \\eta^2 + \\frac{1}{2} h v^2 + \\frac{1}{2} \\hat\\beta v_x^2.
```

`q_global` is a vector of the primitive variables at ALL nodes.
`cache` needs to hold the SBP operators used by the `solver`.

See also [`energy_total_modified`](@ref).
"""
@inline function energy_total_modified!(e, q_global, equations::SvaerdKalischEquations1D,
                                        cache)
    # unpack physical parameters and SBP operator `D1`
    g = gravity(equations)
    (; D1, h, b, v_x, beta_hat) = cache

    # `q_global` is an `ArrayPartition`. It collects the individual arrays for
    # the total water height `eta = h + b` and the velocity `v`.
    eta, v, D = q_global.x
    @.. b = equations.eta0 - D
    @.. h = eta - b

    if D1 isa PeriodicUpwindOperators ||
       D1 isa UpwindOperators
        mul!(v_x, D1.minus, v)
    else
        mul!(v_x, D1, v)
    end

    @.. e = 1 / 2 * g * eta^2 + 1 / 2 * h * v^2 + 1 / 2 * beta_hat * v_x^2
    return e
end
