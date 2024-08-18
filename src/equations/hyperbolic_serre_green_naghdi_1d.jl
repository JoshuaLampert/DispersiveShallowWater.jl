@doc raw"""
    HyperbolicSerreGreenNaghdiEquations1D(bathymetry_type = bathymetry_mild_slope;
                                          gravity_constant, eta0 = 0.0, λ)

Hyperbolic approximation of the Serre-Green-Naghdi system in one spatial
dimension. The equations for flat bathymetry are given by
```math
\begin{aligned}
  h_t + (h v)_x &= 0,\\
  h v_t + \frac{1}{2} g (h^2)_x + \frac{1}{2} h (v^2)_x
    + \biggl( \frac{\lambda}{3} \eta (1 - \eta / h) \biggr)_x &= 0,\\
  h w_t + h v w_x &= \lambda (1 - \eta / h),\\
  \eta_t + \eta_x u &= w.
\end{aligned}
```
The unknown quantities of the hyperbolized Serre-Green-Naghdi equations are the
total water height ``\eta = h + b`` and the velocity ``v``.
The gravitational constant is denoted by `g` and the bottom topography
(bathymetry) ``b = \eta_0 - D``. The water height above the bathymetry
is therefore given by ``h = \eta - \eta_0 + D``.
The total water height is therefore given by ``\eta = h + b``.

Two types of variable `bathymetry_type` are supported:
- [`bathymetry_flat`](@ref): flat bathymetry (typically ``b = 0`` everywhere)
- [`bathymetry_mild_slope`](@ref): variable bathymetry with mild-slope approximation

For the mild-slope approximation, the Serre-Green-Naghdi equations are
```math
\begin{aligned}
  h_t + (h v)_x &= 0,\\
  h v_t + \frac{1}{2} g (h^2)_x + \frac{1}{2} h (v^2)_x
    + \biggl( \frac{\lambda}{3} \eta (1 - \eta / h) \biggr)_x
    + \biggl( g h + \frac{\lambda}{2} (1 - \eta / h) \biggr) b_x &= 0,\\
  h w_t + h v w_x &= \lambda (1 - \eta / h),\\
  \eta_t + \eta_x u + \frac{3}{2} b_x v &= w.
\end{aligned}
```

References for the hyperbolized Serre-Green-Naghdi system can be found in
- Favrie and Gavrilyuk.
  A rapid numerical method for solving Serre-Green-Naghdi equations
  describing long free surface gravity waves
  [DOI: 10.1088/1361-6544/aa712d](https://doi.org/10.1088/1361-6544/aa712d)
- Busto, Dumbser, Escalante, Favrie, and Gavrilyuk.
  On High Order ADER Discontinuous Galerkin Schemes for First Order Hyperbolic
  Reformulations of Nonlinear Dispersive Systems
  [DOI: 10.1007/s10915-021-01429-8](https://doi.org/10.1007/s10915-021-01429-8)

The semidiscretization implemented here conserves
- the total water mass (integral of h) as a linear invariant
- the total modified energy

for periodic boundary conditions (see Ranocha and Ricchiuto (2024)).
Additionally, it is well-balanced for the lake-at-rest stationary solution, see
- Hendrik Ranocha and Mario Ricchiuto (2024)
  Structure-preserving approximations of the Serre-Green-Naghdi
  equations in standard and hyperbolic form
  [arXiv: 2408.02665](https://arxiv.org/abs/2408.02665)
"""
struct HyperbolicSerreGreenNaghdiEquations1D{Bathymetry <: Union{BathymetryFlat, BathymetryMildSlope}, RealT <: Real} <:
       AbstractSerreGreenNaghdiEquations{1, 5}
    bathymetry_type::Bathymetry # type of bathymetry
    gravity::RealT # gravitational constant
    eta0::RealT # constant still-water surface
    λ::RealT # hyperbolic relaxation parameter (→ ∞ for Serre-Green-Naghdi)
end

function HyperbolicSerreGreenNaghdiEquations1D(bathymetry_type = bathymetry_mild_slope;
                                               gravity_constant,
                                               eta0 = 0.0,
                                               λ)
    HyperbolicSerreGreenNaghdiEquations1D(bathymetry_type, gravity_constant, eta0, λ)
end

# TODO: How shall we handle eta=h+b vs the new variable η?
function varnames(::typeof(prim2prim), ::HyperbolicSerreGreenNaghdiEquations1D)
  return ("η", "v", "D", "w", "η")
end
function varnames(::typeof(prim2cons), ::HyperbolicSerreGreenNaghdiEquations1D)
  return ("h", "hv", "b", "hw", "hη")
end

# TODO: There is another name clash. For the SerreGreenNaghdiEquations1D,
#       the corresponding function is called initial_condition_convergence_test
#       However, we cannot use that name since it's not an analytical solution.
#       How shall we handle this?
"""
    initial_condition_soliton(x, t, equations::HyperbolicSerreGreenNaghdiEquations1D, mesh)

A soliton solution of the [`SerreGreenNaghdiEquations1D`](@ref)
used for convergence tests in a periodic domain. This is physically the
same as [`initial_condition_convergence_test`](@ref) for the
[`SerreGreenNaghdiEquations1D`](@ref). Please note that this is not an
exact solution of the [`HyperbolicSerreGreenNaghdiEquations1D`](@ref)
(only in the limit of the relaxation parameter ``\\lambda \\to \\infty``).

See also [`initial_condition_convergence_test`](@ref).
"""
function initial_condition_soliton(x, t, equations::HyperbolicSerreGreenNaghdiEquations1D,
                                            mesh)
    g = gravity_constant(equations)

    # setup parameters data
    h1 = 1.0
    h2 = 1.2
    c = sqrt(g * h2)

    x_t = mod(x - c * t - xmin(mesh), xmax(mesh) - xmin(mesh)) + xmin(mesh)

    h = h1 + (h2 - h1) * sech(x_t / 2 * sqrt(3 * (h2 - h1) / (h1^2 * h2)))^2
    v = c * (1 - h1 / h)
    # w = -h v_x
    w = -h1*sqrt(g*h2)*sqrt((-3*h1 + 3*h2)/(h1^2*h2))*(-h1 + h2)*(-h1 - (-h1 + h2)*sech(x*sqrt((-3*h1 + 3*h2)/(h1^2*h2))/2)^2)*tanh(x*sqrt((-3*h1 + 3*h2)/(h1^2*h2))/2)*sech(x*sqrt((-3*h1 + 3*h2)/(h1^2*h2))/2)^2/(h1 + (-h1 + h2)*sech(x*sqrt((-3*h1 + 3*h2)/(h1^2*h2))/2)^2)^2
    η = h

    return SVector(h, v, 0, w, η)
end

"""
    initial_condition_manufactured(x, t, equations::HyperbolicSerreGreenNaghdiEquations1D, mesh)

A smooth manufactured solution in combination with
[`source_terms_manufactured`](@ref).
"""
function initial_condition_manufactured(x, t, equations::HyperbolicSerreGreenNaghdiEquations1D,
                                            mesh)
    eta = 2 + cospi(2 * (x - 2 * t))
    b = -5 - 2 * cospi(2 * x)
    h = eta - b
    v = sinpi(2 * (x - t / 2))
    D = equations.eta0 - b
    w = 0.0 # TODO: fix, w = -h v_x
    h_approx = h
    return SVector(eta, v, D, w, h_approx)
end

"""
    source_terms_manufactured(q, x, t, equations::HyperbolicSerreGreenNaghdiEquations1D, mesh)

A smooth manufactured solution in combination with
[`initial_condition_manufactured`](@ref).
"""
function source_terms_manufactured(q, x, t, equations::HyperbolicSerreGreenNaghdiEquations1D)
    g = gravity_constant(equations)

    a1 = sinpi(4 * t - 2 * x)
    a2 = cospi(4 * t - 2 * x)
    a3 = sinpi(2 * t - x)
    a4 = cospi(2 * t - x)
    a5 = sinpi(t - 2 * x)
    a6 = cospi(t - 2 * x)
    a7 = sinpi(2 * x)
    a8 = cospi(2 * x)
    a9 = sinpi(2 * t - 4 * x)
    e1 = exp(t)
    e2 = exp(t / 2)

    # Source terms for variable bathymetry
    dh = -4*pi*a1 - a5*(2*pi*a1 - 4*pi*a7) + 2*pi*a6*(a2 + 2*a8 + 7)
    dv = -2*pi*a5*a6 - pi*a6 + 4*pi*a7*g + g*(2*pi*a1 - 4*pi*a7)
    dD = 0.0
    dw = 8*pi^2*a1*a6 - a5*(4*pi^2*a5*(-a2 - 2*a8 - 7) + 2*pi*a6*(-2*pi*a1 + 4*pi*a7)) - 2*pi^2*a5*(-a2 - 2*a8 - 7)
    dh_approx = -4*pi*a1 - 6*pi*a5*a7 - a5*(2*pi*a1 - 4*pi*a7) - 2*pi*a6*(-a2 - 2*a8 - 7)

    return SVector(dh, dv, dD, dw, dh_approx)
end

function create_cache(mesh, equations::HyperbolicSerreGreenNaghdiEquations1D,
                      solver, initial_condition,
                      ::BoundaryConditionPeriodic,
                      RealT, uEltype)
    # TODO: Reduce size of the cache for flat bathymetry
    h = ones(RealT, nnodes(mesh))
    b = zero(h)
    b_x = zero(h)
    η_over_h = zero(h)
    h_x = zero(h)
    v_x = zero(h)
    hv_x = zero(h)
    v2_x = zero(h)
    h_hpb_x = zero(h)
    η_x = zero(h)
    η2_h_x = zero(h)
    w_x = zero(h)
    hvw_x = zero(h)
    tmp = zero(h)

    cache = (; h, b, b_x, η_over_h, h_x, v_x, hv_x, v2_x, h_hpb_x, η_x, η2_h_x, w_x, hvw_x, tmp)
    return cache
end

# Discretization that conserves
# - the total water mass (integral of h) as a linear invariant
# - the total modified energy
# for periodic boundary conditions, see
# - Hendrik Ranocha and Mario Ricchiuto (2024)
#   Structure-preserving approximations of the Serre-Green-Naghdi
#   equations in standard and hyperbolic form
#   [arXiv: 2408.02665](https://arxiv.org/abs/2408.02665)
function rhs!(dq, q, t, mesh,
              equations::HyperbolicSerreGreenNaghdiEquations1D,
              initial_condition,
              ::BoundaryConditionPeriodic,
              source_terms,
              solver, cache)
    # Unpack physical parameters and SBP operator `D1`
    g = gravity_constant(equations)
    (; λ) = equations
    (; D1) = solver

    # `q` and `dq` are `ArrayPartition`s. They collect the individual
    # arrays for the total water height `eta = h + b`, the velocity `v`,
    # and the additional variables `w` and `η`.
    eta, v, D, w, η = q.x
    dh, dv, dD, dw, dη = dq.x # dh = deta since b is constant in time
    fill!(dD, zero(eltype(dD)))

    # TODO: Improve performance for flat bathymetry
    @trixi_timeit timer() "hyperbolic terms" begin
        # Compute all derivatives required below
        (; h, b, b_x, η_over_h, h_x, v_x, hv_x, v2_x, h_hpb_x, η_x, η2_h_x, w_x, hvw_x, tmp) = cache

        @. b = equations.eta0 - D
        @. h = eta - b
        mul!(b_x, D1, b)

        # h_x = D1 * D1
        mul!(h_x, D1, h)

        # v_x = D1 * v
        mul!(v_x, D1, v)

        # hv2_x = D1 * (h * v)
        @. tmp = h * v
        mul!(hv_x, D1, tmp)

        # v2_x = D1 * (v.^2)
        @. tmp = v^2
        mul!(v2_x, D1, tmp)

        # h_hpb_x = D1 * (h .* (h + b)
        @. tmp = h * (h + b)
        mul!(h_hpb_x, D1, tmp)

        # η_x = D1 * η
        mul!(η_x, D1, η)

        # η2_h_x = D1 * (η^2 / h)
        @. η_over_h = η / h
        @. tmp = η * η_over_h
        mul!(η2_h_x, D1, tmp)

        # w_x = D1 * w
        mul!(w_x, D1, w)

        # hvw_x = D1 * (h * v * w)
        @. tmp = h * v * w
        mul!(hvw_x, D1, tmp)

        # Plain: h_t + (h v)_x = 0
        #
        # Split form for energy conservation:
        # h_t + h_x v + h v_x = 0
        @. dh = -(h_x * v + h * v_x)

        # Plain: h v_t + h v v_x + g (h + b) h_x
        #              + ... = 0
        #
        # Split form for energy conservation:
        # h v_t + g (h (h + b))_x - g (h + b) h_x
        #       + 1/2 h (v^2)_x - 1/2 v^2 h_x  + 1/2 v (h v)_x - 1/2 h v v_x
        #       + λ/6 η^2 / h^2 h_x + λ/3 η_x - λ/3 η/h η_x - λ/6 (η^2 / h)_x
        #       + λ/2 b_x - λ/2 η/h b_x = 0
        λ_6 = λ / 6
        λ_3 = λ / 3
        λ_2 = λ / 2
        @. dv = -(g * h_hpb_x - g * (h + b) * h_x
                    + 0.5 * h * v2_x - 0.5 * v^2 * h_x
                    + 0.5 * hv_x * v - 0.5 * h * v * v_x
                    + λ_6 * (η_over_h * η_over_h * h_x - η2_h_x)
                    + λ_3 * (1 - η_over_h) * η_x
                    + λ_2 * (1 - η_over_h) * b_x) / h

        # Plain: h w_t + h v w_x = λ - λ η / h
        #
        # Split form for energy conservation:
        # h w_t + 1/2 (h v w)_x + 1/2 h v w_x
        #       - 1/2 h_x v w - 1/2 h w v_x = λ - λ η / h
        @. dw = ( -(  0.5 * hvw_x
                    + 0.5 * h * v * w_x
                    - 0.5 * h_x * v * w
                    - 0.5 * h * w * v_x) + λ * (1 - η_over_h)) / h

        # No special split form for energy conservation required:
        # η_t + v η_x + 3/2 v b_x = w
        @. dη = -v * η_x - 1.5 * v * b_x + w
    end

    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)

    return nothing
end

@inline function prim2cons(q, equations::HyperbolicSerreGreenNaghdiEquations1D)
    h = waterheight(q, equations)
    v = velocity(q, equations)
    b = bathymetry(q, equations)

    hv = h * v
    hw = h * q[4]
    hη = h * q[5]
    return SVector(h, hv, b, hw, hη)
end

@inline function cons2prim(u, equations::HyperbolicSerreGreenNaghdiEquations1D)
    h, hv, b = u

    eta = h + b
    v = hv / h
    D = equations.eta0 - b
    w = u[4] / h
    η = u[5] / h
    return SVector(eta, v, D, w, η)
end

# The entropy/energy takes the whole `q` for every point in space
"""
    energy_total_modified(q_global, equations::HyperbolicSerreGreenNaghdiEquations1D, cache)

Return the modified total energy of the primitive variables `q_global` for the
[`HyperbolicSerreGreenNaghdiEquations1D`](@ref).
It contains additional terms compared to the usual [`energy_total`](@ref)
modeling non-hydrostatic contributions. The `energy_total_modified`
is a conserved quantity (for periodic boundary conditions).

For a [`bathymetry_mild_slope`](@ref) (and a [`bathymetry_flat`](@ref)),
the total modified energy is given by
```math
\\frac{1}{2} g \\eta^2 + \\frac{1}{2} h v^2 +
\\frac{1}{6} h w^2 + \\frac{\\lambda}{6} h (1 - \\eta / h)^2.
```

`q_global` is a vector of the primitive variables at ALL nodes.
"""
function energy_total_modified(q_global,
                               equations::HyperbolicSerreGreenNaghdiEquations1D,
                               cache)
    # unpack physical parameters and SBP operator `D1`
    g = gravity_constant(equations)
    (; λ) = equations
    (; h, b) = cache

    # `q_global` is an `ArrayPartition`. It collects the individual arrays for
    # the total water height `eta = h + b` and the velocity `v`.
    eta, v, D, w, η = q_global.x
    @. b = equations.eta0 - D
    @. h = eta - b

    e = zero(h)

    # 1/2 g eta^2 + 1/2 h v^2 + 1/6 h^3 w^2 + λ/6 h (1 - η/h)^2
    @. e = 1/2 * g * eta^2 + 1/2 * h * v^2 + 1/6 * h * w^2 + λ/6 * h * (1 - η/h)^2

    return e
end
