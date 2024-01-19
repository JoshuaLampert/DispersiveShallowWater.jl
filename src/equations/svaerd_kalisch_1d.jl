@doc raw"""
    SvärdKalischEquations1D(gravity, eta0 = 1.0, alpha = 0.0, beta = 0.2308939393939394, gamma = 0.04034343434343434)

Dispersive system by Svärd and Kalisch in one spatial dimension with spatially varying bathymetry. The equations are given in conservative variables by
```math
\begin{aligned}
  h_t + (hv)_x &= (\hat\alpha(\hat\alpha(h + b)_x)_x)_x,\\
  (hv)_t + (hv^2)_x + gh(h + b)_x &= (\hat\alpha v(\hat\alpha(h + b)_x)_x)_x + (\hat\beta v_x)_{xt} + \frac{1}{2}(\hat\gamma v_x)_{xx} + \frac{1}{2}(\hat\gamma v_{xx})_x,
\end{aligned}
```
where ``\hat\alpha^2 = \alpha\sqrt{gd}d^2``, ``\hat\beta = \beta d^3``, ``\hat\gamma = \gamma\sqrt{gd}d^3``. The coefficients ``\alpha``, ``\beta`` and ``\gamma`` are provided in dimensionless form and ``d = \eta_0 - b`` is the still-water depth and `eta0` is the still-water surface (lake-at-rest).
The equations can be rewritten in primitive variables as
```math
\begin{aligned}
  \eta_t + ((\eta + D)v)_x = (\hat\alpha(\hat\alpha\eta_x)_x)_x,\\
  v_t(\eta + D) - v((\eta + D)v)_x + ((\eta + D)v^2)_x + g(\eta + D)\eta_x &= (\hat\alpha v(\hat\alpha\eta_x)_x)_x - v(\hat\alpha(\hat\alpha\eta_x)_x)_x + (\hat\beta v_x)_{xt} + \frac{1}{2}(\hat\gamma v_x)_{xx} + \frac{1}{2}(\hat\gamma v_{xx})_x.
\end{aligned}
```
The unknown quantities of the Svärd-Kalisch equations are the total water height ``\eta`` and the velocity ``v``.
The gravitational constant is denoted by `g` and the bottom topography (bathymetry) ``b = -D``. The water height above the bathymetry is therefore given by
``h = \eta + D``.

The equations by Svärd and Kalisch are presented and analyzed in
- Magnus Svärd, Henrik Kalisch (2023)
  A novel energy-bounded Boussinesq model and a well-balanced and stable numerical discretization
  [arXiv: 2302.09924](https://arxiv.org/abs/2302.09924)

"""
struct SvaerdKalischEquations1D{RealT <: Real} <: AbstractSvaerdKalischEquations{1, 3}
    gravity::RealT # gravitational constant
    eta0::RealT    # constant "lake-at-rest" total water height
    alpha::RealT   # coefficient
    beta::RealT    # coefficient
    gamma::RealT   # coefficient
end

const SvärdKalischEquations1D = SvaerdKalischEquations1D

function SvaerdKalischEquations1D(; gravity_constant, eta0 = 1.0, alpha = 0.0,
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
    eta0 = 0.8
    A = 0.02
    # omega = 2*pi/(2.02*sqrt(2))
    k = 0.8406220896381442 # precomputed result of find_zero(k -> omega^2 - equations.gravity * k * tanh(k * eta0), 1.0) using Roots.jl
    if x < -30.5 * pi / k || x > -8.5 * pi / k
        h = 0.0
    else
        h = A * cos(k * x)
    end
    v = sqrt(equations.gravity / k * tanh(k * eta0)) * h / eta0
    if 11.01 <= x && x < 23.04
        b = 0.6 * (x - 11.01) / (23.04 - 11.01)
    elseif 23.04 <= x && x < 27.04
        b = 0.6
    elseif 27.04 <= x && x < 33.07
        b = 0.6 * (33.07 - x) / (33.07 - 27.04)
    else
        b = 0.0
    end
    eta = h + eta0
    D = -b
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
#     D = 3.0
    D = 4 * exp(-(x - 1 / 2)^2)
    return SVector(eta, v, D)
end

"""
    source_terms_manufactured(q, x, t, equations::SvaerdKalischEquations1D, mesh)

A smooth manufactured solution in combination with [`initial_condition_manufactured`](@ref).
"""
function source_terms_manufactured(q, x, t, equations::SvaerdKalischEquations1D)
    g = equations.gravity
#     D = q[3, 1] # D is constant, thus simply take the first entry
    eta0 = equations.eta0
    alpha = equations.alpha
    beta = equations.beta
    gamma = equations.gamma
    a1 = cospi(t - 2 * x)
    a2 = sinpi(t - 2 * x)
    a3 = cospi(4 * t - 2 * x)
    a4 = sinpi(4 * t - 2 * x)

    # Constant bathymetry:
#     dq1 = 8 * pi^3 * alpha * sqrt(g * (D + eta0)) * (D + eta0)^2 * exp(t) * a4 +
#           2 * pi * (D + exp(t) * a3) * exp(t / 2) * a1 - 2 * pi * exp(3 * t / 2) * a2 * a4 -
#           4 * pi * exp(t) * a4 + exp(t) * a3
    # When written in conservative variables:
#     dq2 = 2 * pi * D * g * exp(t) * a4 - D * exp(t / 2) * a2 / 2 -
#           pi * D * exp(t / 2) * a1 - 2 * pi * D * exp(t) * a2 * a1 +
#           8 * pi^3 * alpha * (D + eta0)^2 * sqrt(D * g + eta0 * g) * exp(3 * t / 2) * a1 *
#           a3 - 2 * pi^2 * beta * (D + eta0)^3 * exp(t / 2) * a2 -
#           4 * pi^3 * beta * (D + eta0)^3 * exp(t / 2) * a1 +
#           2 * pi * g * exp(2 * t) * a4 * a3 +
#           8.0 * pi^3 * gamma * (D + eta0)^3 * sqrt(D * g + eta0 * g) * exp(t / 2) * a1 -
#           exp(3 * t / 2) * a2 * a3 / 2 - pi * exp(3 * t / 2) * a1 * a3 -
#           2 * pi * exp(2 * t) * a2 * a1 * a3
#     # When written in primitive variables:
#       dq2 = 8*pi^3*alpha*sqrt(g*(D + eta0))*(D + eta0)^2*exp(3*t/2)*a1*a3 - 2*pi^2*beta*(D + eta0)^3*(a2 + 2*pi*a1)*exp(t/2) + 2*pi*g*(D + exp(t)*a3)*exp(t)*a4 + 8*pi^3*gamma*sqrt(g*(D + eta0))*(D + eta0)^3*exp(t/2)*a1 - (D + exp(t)*a3)*(a2 + 2*pi*a1)*exp(t/2)/2 - 3*pi*(D + exp(t)*a3)*exp(t)*a2*a1 + pi*((D + exp(t)*a3)*a1 - exp(t)*a2*a4)*exp(t)*a2 + pi*exp(2*t)*a2^2*a4

    # Variable bathymetry:
    dq1 = (-10*pi*alpha*sqrt(g*(eta0*exp((x - 0.5)^2) + 4))*(2*x - 1.0)*((10*x - 5.0)*a4 + 2*pi*(eta0*exp((x - 0.5)^2) + 4)*a3)*exp(t + 29*(x - 0.5)^2/4) + 4*pi*alpha*sqrt(g*(eta0*exp((x - 0.5)^2) + 4))*(-10*(x - 0.5)^2*a4 + 2*pi^2*(eta0*exp((x - 0.5)^2) + 4)^2*a4 + 5*(eta0*exp((x - 0.5)^2) + 4)*(-2*(x - 0.5)^2*a4 - 2*pi*(2*x - 1.0)*a3 + a4))*exp(t + 29*(x - 0.5)^2/4) + (-4*pi*a4 + a3)*exp(t + 39*(x - 0.5)^2/4) + (2*pi*(exp(t + (x - 0.5)^2)*a3 + 4)*a1 + (8*x - 2*pi*exp(t + (x - 0.5)^2)*a4 - 4.0)*a2)*exp(t/2 + 35*(x - 0.5)^2/4))*exp(-39*(x - 0.5)^2/4)

    # When written in conservative variables:
#         dq2 = (-(10 * pi * alpha * sqrt(g * (eta0 * exp((x - 0.5)^2) + 4)) * (2 * x - 1.0) *
#                  ((10 * x - 5.0) * a4 + 2 * pi * (eta0 * exp((x - 0.5)^2) + 4) * a3) *
#                  exp(t + 29 * (x - 0.5)^2 / 4) +
#                  4 * pi * alpha * sqrt(g * (eta0 * exp((x - 0.5)^2) + 4)) *
#                  (10 * (x - 0.5)^2 * a4 - 2 * pi^2 * (eta0 * exp((x - 0.5)^2) + 4)^2 * a4 +
#                   5 * (eta0 * exp((x - 0.5)^2) + 4) *
#                   (2 * (x - 0.5)^2 * a4 + 2 * pi * (2 * x - 1.0) * a3 - a4)) *
#                  exp(t + 29 * (x - 0.5)^2 / 4) +
#                  (4 * pi * a4 - a3) * exp(t + 39 * (x - 0.5)^2 / 4) +
#                  (-2 * pi * (exp(t + (x - 0.5)^2) * a3 + 4) * a1 +
#                   (-8 * x + 2 * pi * exp(t + (x - 0.5)^2) * a4 + 4.0) * a2) *
#                  exp(t / 2 + 35 * (x - 0.5)^2 / 4)) * exp(t / 2 + 20 * (x - 0.5)^2) * a2 +
#                (20 * pi * alpha * sqrt(g * (eta0 * exp((x - 0.5)^2) + 4)) * (2 * x - 1.0) *
#                 ((10 * x - 5.0) * a4 + 2 * pi * (eta0 * exp((x - 0.5)^2) + 4) * a3) *
#                 exp(3 * t / 2 + 35 * (x - 0.5)^2 / 2) * a2 +
#                 8 * pi^2 * alpha * sqrt(g * (eta0 * exp((x - 0.5)^2) + 4)) *
#                 (eta0 * exp((x - 0.5)^2) + 4) *
#                 ((10 * x - 5.0) * a4 + 2 * pi * (eta0 * exp((x - 0.5)^2) + 4) * a3) *
#                 exp(3 * t / 2 + 35 * (x - 0.5)^2 / 2) * a1 +
#                 8 * pi * alpha * sqrt(g * (eta0 * exp((x - 0.5)^2) + 4)) *
#                 (10 * (x - 0.5)^2 * a4 - 2 * pi^2 * (eta0 * exp((x - 0.5)^2) + 4)^2 * a4 +
#                  5 * (eta0 * exp((x - 0.5)^2) + 4) *
#                  (2 * (x - 0.5)^2 * a4 + 2 * pi * (2 * x - 1.0) * a3 - a4)) *
#                 exp(3 * t / 2 + 35 * (x - 0.5)^2 / 2) * a2 -
#                 4 * pi * beta * (eta0 * exp((x - 0.5)^2) + 4)^2 *
#                 ((6.0 - 12 * x) * a1 + 12 * pi * (2 * x - 1.0) * a2 +
#                  pi * (eta0 * exp((x - 0.5)^2) + 4) * a2 +
#                  2 * pi^2 * (eta0 * exp((x - 0.5)^2) + 4) * a1) *
#                 exp(t / 2 + 17 * (x - 0.5)^2) +
#                 pi * gamma * sqrt(g * (eta0 * exp((x - 0.5)^2) + 4)) *
#                 (eta0 * exp((x - 0.5)^2) + 4) *
#                 (-384.0 * (x - 0.5)^2 * a1 +
#                  168.0 * pi * (2 * x - 1.0) * (eta0 * exp((x - 0.5)^2) + 4) * a2 +
#                  16.0 * pi^2 * (eta0 * exp((x - 0.5)^2) + 4)^2 * a1 +
#                  (32.0 * (x - 0.5)^2 -
#                   8.0 * (eta0 * exp((x - 0.5)^2) + 4) * (2 * (x - 0.5)^2 - 1)) * a1 +
#                  (48.0 * eta0 * exp((x - 0.5)^2) -
#                   (x - 0.5)^2 * (96.0 * eta0 * exp((x - 0.5)^2) + 384.0) - 768.0 * (x - 0.5)^2 +
#                   192.0) * a1) * exp(t / 2 + 33 * (x - 0.5)^2 / 2) +
#                 2 * (4 * pi * a4 - a3) * exp(3 * t / 2 + 20 * (x - 0.5)^2) * a2 +
#                 (4 * pi * g * (exp(t + (x - 0.5)^2) * a3 + 4) * exp(t) * a4 -
#                  (exp(t + (x - 0.5)^2) * a3 + 4) * exp(t / 2) * a2 -
#                  2 * pi * (exp(t + (x - 0.5)^2) * a3 + 4) * exp(t / 2) * a1 -
#                  8 * pi * (exp(t + (x - 0.5)^2) * a3 + 4) * exp(t) * a2 * a1 +
#                  (-16 * x + 4 * pi * exp(t + (x - 0.5)^2) * a4 + 8.0) * exp(t) * a2^2) *
#                 exp(19 * (x - 0.5)^2)) * exp(39 * (x - 0.5)^2 / 4) / 2) *
#               exp(-119 * (x - 0.5)^2 / 4)# - 180.0
    # When written in primitive variables:
    dq2 = (5*pi*alpha*sqrt(g*(eta0*exp((x - 0.5)^2) + 4))*(2*x - 1.0)*((10*x - 5.0)*a4 + 2*pi*(eta0*exp((x - 0.5)^2) + 4)*a3)*exp(3*t/2 + 105*(x - 0.5)^2/4)*a2 + 4*pi^2*alpha*sqrt(g*(eta0*exp((x - 0.5)^2) + 4))*(eta0*exp((x - 0.5)^2) + 4)*((10*x - 5.0)*a4 + 2*pi*(eta0*exp((x - 0.5)^2) + 4)*a3)*exp(3*t/2 + 105*(x - 0.5)^2/4)*a1 + pi*alpha*sqrt(g*(eta0*exp((x - 0.5)^2) + 4))*(-5*(2*x - 1.0)*((10*x - 5.0)*a4 + 2*pi*(eta0*exp((x - 0.5)^2) + 4)*a3)*exp(t + 25*(x - 0.5)^2/4) + 2*(-10*(x - 0.5)^2*a4 + 2*pi^2*(eta0*exp((x - 0.5)^2) + 4)^2*a4 + 5*(eta0*exp((x - 0.5)^2) + 4)*(-2*(x - 0.5)^2*a4 - 2*pi*(2*x - 1.0)*a3 + a4))*exp(t + 25*(x - 0.5)^2/4))*exp(t/2 + 20*(x - 0.5)^2)*a2 + 2*pi*alpha*sqrt(g*(eta0*exp((x - 0.5)^2) + 4))*(10*(x - 0.5)^2*a4 - 2*pi^2*(eta0*exp((x - 0.5)^2) + 4)^2*a4 + 5*(eta0*exp((x - 0.5)^2) + 4)*(2*(x - 0.5)^2*a4 + 2*pi*(2*x - 1.0)*a3 - a4))*exp(3*t/2 + 105*(x - 0.5)^2/4)*a2 + 2*pi*beta*(eta0*exp((x - 0.5)^2) + 4)^2*((6.0 - 12*x)*(2*pi*a2 - a1) - pi*(eta0*exp((x - 0.5)^2) + 4)*(a2 + 2*pi*a1))*exp(t/2 + 103*(x - 0.5)^2/4) + 4*pi*gamma*sqrt(g*(eta0*exp((x - 0.5)^2) + 4))*(eta0*exp((x - 0.5)^2) + 4)*(-48*(x - 0.5)^2*a1 + 21*pi*(2*x - 1.0)*(eta0*exp((x - 0.5)^2) + 4)*a2 + 2*pi^2*(eta0*exp((x - 0.5)^2) + 4)^2*a1 + ((1 - 2*(x - 0.5)^2)*(eta0*exp((x - 0.5)^2) + 4) + 4*(x - 0.5)^2)*a1 + 6*(eta0*exp((x - 0.5)^2) - 2*(x - 0.5)^2*(eta0*exp((x - 0.5)^2) + 4) - 16*(x - 0.5)^2 + 4)*a1)*exp(t/2 + 101*(x - 0.5)^2/4) + (4*pi*g*(exp(t + (x - 0.5)^2)*a3 + 4)*exp(t)*a4 - (exp(t + (x - 0.5)^2)*a3 + 4)*(a2 + 2*pi*a1)*exp(t/2) - 6*pi*(exp(t + (x - 0.5)^2)*a3 + 4)*exp(t)*a2*a1 + (2*pi*(exp(t + (x - 0.5)^2)*a3 + 4)*a1 + (8*x - 2*pi*exp(t + (x - 0.5)^2)*a4 - 4.0)*a2)*exp(t)*a2 + (-8*x + 2*pi*exp(t + (x - 0.5)^2)*a4 + 4.0)*exp(t)*a2^2)*exp(111*(x - 0.5)^2/4)/2)*exp(-115*(x - 0.5)^2/4)

    return SVector(dq1, dq2, 0.0)
end

function create_cache(mesh,
                      equations::SvaerdKalischEquations1D,
                      solver,
                      initial_condition,
                      RealT,
                      uEltype)
    #  Assume D is independent of time and compute D evaluated at mesh points once.
    D = Array{RealT}(undef, nnodes(mesh))
    x = grid(solver)
    for i in eachnode(solver)
        D[i] = initial_condition(x[i], 0.0, equations, mesh)[3]
    end
    d = equations.eta0 .+ D
    h = Array{RealT}(undef, nnodes(mesh))
    hv = similar(h)
    alpha_hat = sqrt.(equations.alpha * sqrt.(equations.gravity * d) .* d .^ 2)
    beta_hat = equations.beta * d .^ 3
    gamma_hat = equations.gamma * sqrt.(equations.gravity * d) .* d .^ 3
    tmp1 = similar(h)
    tmp2 = similar(h)
    hmD1betaD1 = Array{RealT}(undef, nnodes(mesh), nnodes(mesh))
    if solver.D1 isa PeriodicDerivativeOperator ||
       solver.D1 isa UniformPeriodicCoupledOperator
        D1_central = solver.D1
        sparse_D1 = sparse(D1_central)
        D1betaD1 = sparse_D1 * Diagonal(beta_hat) * sparse_D1
    elseif solver.D1 isa PeriodicUpwindOperators
        D1_central = solver.D1.central
        D1betaD1 = sparse(solver.D1.plus) * Diagonal(beta_hat) * sparse(solver.D1.minus)
    else
        @error "unknown type of first-derivative operator"
    end
    return (hmD1betaD1 = hmD1betaD1, D1betaD1 = D1betaD1, d = d, h = h, hv = hv,
            alpha_hat = alpha_hat, beta_hat = beta_hat, gamma_hat = gamma_hat,
            tmp1 = tmp1, tmp2 = tmp2, D1_central = D1_central, D1 = solver.D1)
end

# Discretization that conserves the mass (for eta and v) and is energy-bounded for periodic boundary conditions
function rhs!(du_ode, u_ode, t, mesh, equations::SvaerdKalischEquations1D,
              initial_condition, ::BoundaryConditionPeriodic, source_terms,
              solver, cache)
    @unpack hmD1betaD1, D1betaD1, d, h, hv, alpha_hat, beta_hat, gamma_hat, tmp1, tmp2, D1_central = cache
    q = wrap_array(u_ode, mesh, equations, solver)
    dq = wrap_array(du_ode, mesh, equations, solver)

    eta = view(q, 1, :)
    v = view(q, 2, :)
    D = view(q, 3, :)
    deta = view(dq, 1, :)
    dv = view(dq, 2, :)
    dD = view(dq, 3, :)
    fill!(dD, zero(eltype(dD)))

    @timeit timer() "deta hyperbolic" begin
        h = eta .+ D
        hv = h .* v

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
            @error "unknown type of first derivative operator"
        end
    end

    # split form
    @timeit timer() "dv hyperbolic" begin
        dv[:] = -(0.5 * (D1_central * (hv .* v) + hv .* D1v - v .* (D1_central * hv)) +
                  equations.gravity * h .* D1eta +
                  0.5 * (vD1y - D1vy - yD1v) -
                  0.5 * D1_central * (gamma_hat .* (solver.D2 * v)) -
                  0.5 * solver.D2 * (gamma_hat .* D1v))
    end

    # no split form
    #     dv[:] = -(D1_central * (hv .* v) - v .* (D1_central * hv)+
    #               equations.gravity * h .* D1eta +
    #               vD1y - D1vy -
    #               0.5 * D1_central * (gamma_hat .* (solver.D2 * v)) -
    #               0.5 * solver.D2 * (gamma_hat .* D1v))

    @timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations, solver)
    @timeit timer() "dv elliptic" begin
        hmD1betaD1 = Diagonal(h) - D1betaD1
        dv[:] = hmD1betaD1 \ dv
    end

    return nothing
end

@inline function prim2cons(q, equations::SvaerdKalischEquations1D)
    eta, v, D = q

    h = eta + D
    hv = h * v
    b = -D
    return SVector(h, hv, b)
end

@inline function cons2prim(u, equations::SvaerdKalischEquations1D)
    h, hv, b = u

    eta = h + b
    v = hv / h
    D = -b
    return SVector(eta, v, D)
end

@inline function waterheight_total(q, equations::SvaerdKalischEquations1D)
    return q[1]
end

@inline function velocity(q, equations::SvaerdKalischEquations1D)
    return q[2]
end

@inline function bathymetry(q, equations::SvaerdKalischEquations1D)
    return -q[3]
end

@inline function waterheight(q, equations::SvaerdKalischEquations1D)
    return waterheight_total(q, equations) - bathymetry(q, equations)
end

@inline function energy_total(q, equations::SvaerdKalischEquations1D)
    eta, v, D = q
    e = 0.5 * (equations.gravity * eta^2 + (D + eta) * v^2)
    return e
end

@inline entropy(u, equations::SvaerdKalischEquations1D) = energy_total(u, equations)

# The modified entropy/total energy takes the whole `q` for every point in space
"""
    energy_total_modified(q, equations::SvaerdKalischEquations1D, cache)

Return the modified total energy of the primitive variables `q` for the
`SvaerdKalischEquations1D`. It contains an additional term containing a
derivative compared to the usual `energy_total`. The `energy_total_modified`
is a conserved quantity of the Svärd-Kalisch equations.

`q` is a vector of the primitive variables at ALL nodes, i.e., a matrix
of the correct length `nvariables(equations)` as first dimension and the
number of nodes as length of the second dimension.
`cache` needs to hold the first-derivative SBP operator `D1`.
"""
@inline function energy_total_modified(q, equations::SvaerdKalischEquations1D, cache)
    e_modified = zeros(eltype(q), size(q, 2))
    # Need to compute new beta_hat, do not use the old one from the `cache`
    eta = view(q, 1, :)
    v = view(q, 2, :)
    D = view(q, 3, :)
    beta_hat = equations.beta * (eta .+ D) .^ 3
    if cache.D1 isa PeriodicDerivativeOperator ||
       cache.D1 isa UniformPeriodicCoupledOperator
        tmp = 0.5 * beta_hat .* ((cache.D1 * v) .^ 2)
    elseif cache.D1 isa PeriodicUpwindOperators
        tmp = 0.5 * beta_hat .* ((cache.D1.minus * v) .^ 2)
    else
        @error "unknown type of first-derivative operator"
    end
    for i in 1:size(q, 2)
        e_modified[i] = energy_total(view(q, :, i), equations) + tmp[i]
    end
    return e_modified
end

varnames(::typeof(energy_total_modified), equations) = ("e_modified",)

"""
    entropy_modified(q, equations::SvaerdKalischEquations1D, cache)

Alias for [`energy_total_modified`](@ref).
"""
@inline function entropy_modified(q, equations::SvaerdKalischEquations1D, cache)
    energy_total_modified(q, equations, cache)
end

varnames(::typeof(entropy_modified), equations) = ("U_modified",)

# Calculate the error for the "lake-at-rest" test case where eta should
# be a constant value over time
@inline function lake_at_rest_error(u, equations::SvaerdKalischEquations1D)
    eta, _, _ = u
    return abs(equations.eta0 - eta)
end
