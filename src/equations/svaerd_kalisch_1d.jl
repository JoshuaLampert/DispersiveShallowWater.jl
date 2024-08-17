@doc raw"""
    SvaerdKalischEquations1D(; gravity_constant, eta0 = 0.0,
                               alpha = 0.0,
                               beta = 0.2308939393939394,
                               gamma = 0.04034343434343434)

Dispersive system by Svärd and Kalisch in one spatial dimension with spatially varying bathymetry. The equations are given in conservative variables by
```math
\begin{aligned}
  h_t + (hv)_x &= (\hat\alpha(\hat\alpha(h + b)_x)_x)_x,\\
  (hv)_t + (hv^2)_x + gh(h + b)_x &= (\hat\alpha v(\hat\alpha(h + b)_x)_x)_x + (\hat\beta v_x)_{xt} + \frac{1}{2}(\hat\gamma v_x)_{xx} + \frac{1}{2}(\hat\gamma v_{xx})_x,
\end{aligned}
```
where ``\hat\alpha^2 = \alpha\sqrt{gD}D^2``, ``\hat\beta = \beta D^3``, ``\hat\gamma = \gamma\sqrt{gD}D^3``. The coefficients ``\alpha``, ``\beta`` and ``\gamma`` are provided in dimensionless form and ``D = \eta_0 - b`` is the still-water depth and `eta0` is the still-water surface (lake-at-rest).
The equations can be rewritten in primitive variables as
```math
\begin{aligned}
  \eta_t + ((\eta + D)v)_x &= (\hat\alpha(\hat\alpha\eta_x)_x)_x,\\
  v_t(\eta + D) - v((\eta + D)v)_x + ((\eta + D)v^2)_x + g(\eta + D)\eta_x &= (\hat\alpha v(\hat\alpha\eta_x)_x)_x - v(\hat\alpha(\hat\alpha\eta_x)_x)_x + (\hat\beta v_x)_{xt} + \frac{1}{2}(\hat\gamma v_x)_{xx} + \frac{1}{2}(\hat\gamma v_{xx})_x.
\end{aligned}
```
The unknown quantities of the Svärd-Kalisch equations are the total water height ``\eta`` and the velocity ``v``.
The gravitational constant is denoted by `g` and the bottom topography (bathymetry) ``b = \eta_0 - D``. The water height above the bathymetry is therefore given by
``h = \eta - \eta_0 + D``.

`SvärdKalischEquations1D` is an alias for `SvaerdKalischEquations1D`.

The equations by Svärd and Kalisch are presented and analyzed in Svärd and Kalisch (2023).
The semidiscretization implemented here conserves the mass and the energy, is well-balanced for the lake-at-rest state,
and is developed in Lampert and Ranocha (2024).

- Magnus Svärd, Henrik Kalisch (2023)
  A novel energy-bounded Boussinesq model and a well-balanced and stable numerical discretization
  [arXiv: 2302.09924](https://arxiv.org/abs/2302.09924)
- Joshua Lampert, Hendrik Ranocha (2024)
  Structure-Preserving Numerical Methods for Two Nonlinear Systems of Dispersive Wave Equations
  [DOI: 10.48550/arXiv.2402.16669](https://doi.org/10.48550/arXiv.2402.16669)
"""
struct SvaerdKalischEquations1D{RealT <: Real} <: AbstractSvaerdKalischEquations{1, 3}
    gravity::RealT # gravitational constant
    eta0::RealT    # constant still-water surface
    alpha::RealT   # coefficient
    beta::RealT    # coefficient
    gamma::RealT   # coefficient
end

const SvärdKalischEquations1D = SvaerdKalischEquations1D

function SvaerdKalischEquations1D(; gravity_constant, eta0 = 0.0, alpha = 0.0,
                                  beta = 0.2308939393939394, gamma = 0.04034343434343434)
    SvaerdKalischEquations1D(gravity_constant, eta0, alpha, beta, gamma)
end

varnames(::typeof(prim2prim), ::SvaerdKalischEquations1D) = ("η", "v", "D")
varnames(::typeof(prim2cons), ::SvaerdKalischEquations1D) = ("h", "hv", "b")

"""
    initial_condition_dingemans(x, t, equations::SvaerdKalischEquations1D, mesh)

The initial condition that uses the dispersion relation of the Euler equations
to approximate waves generated by a wave maker as it is done by experiments of
Dingemans. The topography is a trapezoidal. It is assumed that `equations.eta0 = 0.8`.

References:
- Magnus Svärd, Henrik Kalisch (2023)
  A novel energy-bounded Boussinesq model and a well-balanced and stable numerical discretization
  [arXiv: 2302.09924](https://arxiv.org/abs/2302.09924)
- Maarten W. Dingemans (1994)
  Comparison of computations with Boussinesq-like models and laboratory measurements
  [link](https://repository.tudelft.nl/islandora/object/uuid:c2091d53-f455-48af-a84b-ac86680455e9/datastream/OBJ/download)
"""
function initial_condition_dingemans(x, t, equations::SvaerdKalischEquations1D, mesh)
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
    eta = h + h0
    D = equations.eta0 - b
    return SVector(eta, v, D)
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
    g = equations.gravity
    eta0 = equations.eta0
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

function create_cache(mesh, equations::SvaerdKalischEquations1D,
                      solver, initial_condition,
                      ::BoundaryConditionPeriodic,
                      RealT, uEltype)
    #  Assume D is independent of time and compute D evaluated at mesh points once.
    D = Array{RealT}(undef, nnodes(mesh))
    x = grid(solver)
    for i in eachnode(solver)
        D[i] = still_waterdepth(initial_condition(x[i], 0.0, equations, mesh), equations)
    end
    h = Array{RealT}(undef, nnodes(mesh))
    hv = similar(h)
    alpha_hat = sqrt.(equations.alpha * sqrt.(equations.gravity * D) .* D .^ 2)
    beta_hat = equations.beta * D .^ 3
    gamma_hat = equations.gamma * sqrt.(equations.gravity * D) .* D .^ 3
    tmp2 = similar(h)
    M = mass_matrix(solver.D1)
    if solver.D1 isa PeriodicDerivativeOperator ||
       solver.D1 isa UniformPeriodicCoupledOperator
        D1_central = solver.D1
        sparse_D1 = sparse(D1_central)
        # We use the periodic SBP property
        #   M * sparse_D1 == -sparse_D1' * M
        # to avoid possible floating point errors in the symmetry of the matrix.
        minus_MD1betaD1 = sparse_D1' * M * Diagonal(beta_hat) * sparse_D1
    elseif solver.D1 isa PeriodicUpwindOperators
        D1_central = solver.D1.central
        # We use the periodic upwind SBP property
        #   M * sparse_D1plus == -sparse_D1minus' * M
        # to avoid possible floating point errors in the symmetry of the matrix.
        sparse_D1minus = sparse(solver.D1.minus)
        minus_MD1betaD1 = sparse_D1minus' * M * Diagonal(beta_hat) * sparse_D1minus
    else
        @error "unknown type of first-derivative operator: $(typeof(solver.D1))"
    end
    factorization = cholesky(Symmetric(M * Diagonal(ones(nnodes(mesh))) + minus_MD1betaD1))
    return (factorization = factorization, minus_MD1betaD1 = minus_MD1betaD1, D = D, h = h,
            hv = hv, alpha_hat = alpha_hat, gamma_hat = gamma_hat,
            tmp2 = tmp2, D1_central = D1_central, M = M, D1 = solver.D1)
end

# Discretization that conserves the mass (for eta and for flat bottom hv) and the energy for periodic boundary conditions, see
# - Joshua Lampert and Hendrik Ranocha (2024)
#   Structure-Preserving Numerical Methods for Two Nonlinear Systems of Dispersive Wave Equations
#   [DOI: 10.48550/arXiv.2402.16669](https://doi.org/10.48550/arXiv.2402.16669)
function rhs!(dq, q, t, mesh, equations::SvaerdKalischEquations1D,
              initial_condition, ::BoundaryConditionPeriodic, source_terms,
              solver, cache)
    @unpack factorization, minus_MD1betaD1, D, h, hv, alpha_hat, gamma_hat, tmp1, tmp2, D1_central, M = cache

    eta = q.x[1]
    v = q.x[2]
    deta = dq.x[1]
    dv = dq.x[2]
    dD = dq.x[3]
    fill!(dD, zero(eltype(dD)))

    @trixi_timeit timer() "deta hyperbolic" begin
        @. h = eta + D - equations.eta0
        @. hv = h * v

        if solver.D1 isa PeriodicDerivativeOperator ||
           solver.D1 isa UniformPeriodicCoupledOperator
            D1eta = D1_central * eta
            D1v = D1_central * v
            tmp1 = alpha_hat .* (D1_central * (alpha_hat .* D1eta))
            vD1y = v .* (D1_central * tmp1)
            D1vy = D1_central * (v .* tmp1)
            yD1v = tmp1 .* D1v
            @. tmp2 = tmp1 - hv
            mul!(deta, D1_central, tmp2)
        elseif solver.D1 isa PeriodicUpwindOperators
            D1eta = D1_central * eta
            D1v = D1_central * v
            tmp1 = alpha_hat .* (solver.D1.minus * (alpha_hat .* (solver.D1.plus * eta)))
            vD1y = v .* (solver.D1.minus * tmp1)
            D1vy = solver.D1.minus * (v .* tmp1)
            yD1v = tmp1 .* (solver.D1.plus * v)
            deta[:] = solver.D1.minus * tmp1 - D1_central * hv
        else
            @error "unknown type of first derivative operator: $(typeof(solver.D1))"
        end
    end

    # split form
    @trixi_timeit timer() "dv hyperbolic" begin
        D1_hv = D1_central * hv
        D1_hv2 = D1_central * (hv .* v)
        D1_gamma_hat_D2_v = D1_central * (gamma_hat .* (solver.D2 * v))
        D2_gamma_hat_D1_v = solver.D2 * (gamma_hat .* D1v)
        @. dv = -(0.5 * (D1_hv2 + hv * D1v - v * D1_hv) +
                  equations.gravity * h * D1eta +
                  0.5 * (vD1y - D1vy - yD1v) -
                  0.5 * D1_gamma_hat_D2_v -
                  0.5 * D2_gamma_hat_D1_v)
    end

    # no split form
    #     dv[:] = -(D1_central * (hv .* v) - v .* (D1_central * hv)+
    #               equations.gravity * h .* D1eta +
    #               vD1y - D1vy -
    #               0.5 * D1_central * (gamma_hat .* (solver.D2 * v)) -
    #               0.5 * solver.D2 * (gamma_hat .* D1v))

    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)
    @trixi_timeit timer() "dv elliptic" begin
        # decompose M * (h - D1betaD1) because it is guaranteed to be symmetric and pos. def.,
        # while (h - D1betaD1) is not necessarily
        hmD1betaD1 = Symmetric(M * Diagonal(h) + minus_MD1betaD1)
        # If the time integration method takes a too large step, h become become negative
        # letting the factorization fail. In this case the solution is set to `NaN` to force
        # rejecting the step
        cholesky!(factorization, hmD1betaD1, check = false)
        if issuccess(factorization)
            mul!(tmp1, M, dv)
            dv[:] = factorization \ tmp1
        else
            fill!(dv, NaN)
        end
    end

    return nothing
end

@inline function prim2cons(q, equations::SvaerdKalischEquations1D)
    eta, v, D = q

    b = bathymetry(q, equations)
    h = eta - b
    hv = h * v
    return SVector(h, hv, b)
end

@inline function cons2prim(u, equations::SvaerdKalischEquations1D)
    h, hv, b = u

    eta = h + b
    v = hv / h
    D = equations.eta0 - b
    return SVector(eta, v, D)
end

@inline function waterheight_total(q, equations::SvaerdKalischEquations1D)
    return q[1]
end

@inline function velocity(q, equations::SvaerdKalischEquations1D)
    return q[2]
end

@inline function bathymetry(q, equations::SvaerdKalischEquations1D)
    D = q[3]
    return equations.eta0 - D
end

@inline function waterheight(q, equations::SvaerdKalischEquations1D)
    return waterheight_total(q, equations) - bathymetry(q, equations)
end

@inline entropy(u, equations::SvaerdKalischEquations1D) = energy_total(u, equations)

# The modified entropy/total energy takes the whole `q` for every point in space
"""
    energy_total_modified(q_global, equations::SvaerdKalischEquations1D, cache)

Return the modified total energy of the primitive variables `q_global` for the
[`SvaerdKalischEquations1D`](@ref). It contains an additional term containing a
derivative compared to the usual [`energy_total`](@ref). The `energy_total_modified`
is a conserved quantity of the Svärd-Kalisch equations given by
```math
\\frac{1}{2} g h^2 + \\frac{1}{2} h v^2 + \\frac{1}{2} \\hat\\beta v_x^2.
```

`q_global` is a vector of the primitive variables at ALL nodes.
`cache` needs to hold the first-derivative SBP operator `D1`.
"""
@inline function energy_total_modified(q_global, equations::SvaerdKalischEquations1D, cache)
    # Need to compute new beta_hat, do not use the old one from the `cache`
    v = q_global.x[2]
    D = q_global.x[3]
    N = length(v)
    e_modified = zeros(eltype(q_global), N)
    beta_hat = equations.beta * D .^ 3
    if cache.D1 isa PeriodicDerivativeOperator ||
       cache.D1 isa UniformPeriodicCoupledOperator
        tmp = 0.5 * beta_hat .* ((cache.D1 * v) .^ 2)
    elseif cache.D1 isa PeriodicUpwindOperators
        tmp = 0.5 * beta_hat .* ((cache.D1.minus * v) .^ 2)
    else
        @error "unknown type of first-derivative operator: $(typeof(cache.D1))"
    end
    for i in 1:N
        e_modified[i] = energy_total(get_node_vars(q_global, equations, i), equations) +
                        tmp[i]
    end
    return e_modified
end

# Calculate the error for the "lake-at-rest" test case where eta should
# be a constant value over time
@inline function lake_at_rest_error(u, equations::SvaerdKalischEquations1D)
    eta, _, _ = u
    return abs(equations.eta0 - eta)
end
