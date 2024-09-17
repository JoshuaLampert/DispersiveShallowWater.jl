@doc raw"""
    BBMEquation1D(; bathymetry_type = bathymetry_flat,
                  gravity_constant, D = 1.0, eta0 = 0.0, split_form = true)

BBM (Benjamin–Bona–Mahony) equation in one spatial dimension. The equations are given by
```math
\begin{aligned}
  \eta_t + \sqrt{gD}\eta_x + \frac{3}{2}\sqrt{\frac{g}{D}}\eta\eta_x - \frac{1}{6}D^2\eta_{xxt} &= 0.
\end{aligned}
```

The unknown quantity of the BBM equation is the total water height ``\eta``.
The gravitational constant is denoted by `g` and the constant bottom topography (bathymetry) ``b = \eta_0 - D``,
where ``\eta_0`` is the constant still-water surface and ``D`` the still-water depth. The water height above
the bathymetry is therefore given by ``h = \eta - \eta_0 + D``. The BBM equation is only implemented for ``\eta_0 = 0``.

Currently, the equations only support a flat bathymetry, see [`bathymetry_flat`](@ref).

The BBM equation is first described in Benjamin, Bona, and Mahony (1972).
The semidiscretization implemented here is developed in Ranocha, Mitsotakis, and Ketcheson (2020)
for `split_form = true` and in Linders, Ranocha, and Birken (2023) for `split_form = false`.
If `split_form` is `true` a split form in the semidiscretization is used, which conserves
- the total water mass (integral of ``h``) as a linear invariant
- a quadratic invariant (integral of ``\eta^2 + \eta_x^2``), which is called here [`energy_total_modified`](@ref)
  (and [`entropy_modified`](@ref)) because it contains derivatives of the solution

for periodic boundary conditions. If `split_form` is `false` the semidiscretization conserves
- the total water mass (integral of ``h``) as a linear invariant
- the Hamiltonian (integral of ``1/6 \eta^3 + 1/2 \eta^2``) (see [`hamiltonian`](@ref))

for periodic boundary conditions.

- Thomas B. Benjamin, Jerry L. Bona and John J. Mahony (1972)
  Model equations for long waves in nonlinear dispersive systems
  [DOI: 10.1098/rsta.1972.0032](https://doi.org/10.1098/rsta.1972.0032)
- Hendrik Ranocha, Dimitrios Mitsotakis and David I. Ketcheson (2020)
  A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations
  [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)
- Viktor Linders, Hendrik Ranocha and Philipp Birken (2023)
  Resolving entropy growth from iterative methods
  [DOI: 10.1007/s10543-023-00992-w](https://doi.org/10.1007/s10543-023-00992-w)
"""
struct BBMEquation1D{Bathymetry <: AbstractBathymetry, RealT <: Real} <:
       AbstractBBMEquation{1, 1}
    bathymetry_type::Bathymetry # type of bathymetry
    gravity::RealT # gravitational constant
    D::RealT # still-water depth
    eta0::RealT # constant still-water surface
    split_form::Bool # whether to use a split-form or not
end

function BBMEquation1D(; bathymetry_type = bathymetry_flat,
                       gravity_constant, D = 1.0, eta0 = 0.0, split_form = true)
    BBMEquation1D(bathymetry_type, gravity_constant, D, eta0, split_form)
end

"""
    initial_condition_convergence_test(x, t, equations::BBMEquation1D, mesh)

A travelling-wave solution used for convergence tests in a periodic domain, here generalized
for dimensional variables.

See section 4.1.3 in (there is an error in paper: it should be `sech^2` instead of `cosh`):
- Hendrik Ranocha, Dimitrios Mitsotakis and David I. Ketcheson (2020)
  A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations
  [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)
"""
function initial_condition_convergence_test(x, t, equations::BBMEquation1D, mesh)
    g = gravity_constant(equations)
    D = equations.D
    alpha = sqrt(g * D)
    beta = 3 / 2 * sqrt(g / D)
    gamma = 1 / 6 * D^2
    c = 1.2 * alpha
    A = 3 * (c - alpha) / beta
    K = 0.5 * sqrt(1 / gamma * (1 - alpha / c))
    x_t = mod(x - c * t - xmin(mesh), xmax(mesh) - xmin(mesh)) + xmin(mesh)
    eta = A * sech(K * x_t)^2
    return SVector(eta)
end

"""
    initial_condition_manufactured(x, t, equations::BBMEquation1D, mesh)

A smooth manufactured solution in combination with [`source_terms_manufactured`](@ref).
"""
function initial_condition_manufactured(x, t,
                                        equations::BBMEquation1D,
                                        mesh)
    eta = exp(t / 2) * sinpi(2 * (x - t / 2))
    return SVector(eta)
end

"""
    source_terms_manufactured(q, x, t, equations::BBMEquation1D, mesh)

A smooth manufactured solution in combination with [`initial_condition_manufactured`](@ref).
"""
function source_terms_manufactured(q, x, t, equations::BBMEquation1D)
    g = gravity_constant(equations)
    D = still_waterdepth(q, equations)
    a1 = sqrt(g * D)
    a2 = sqrt(g / D)
    a3 = cospi(t - 2 * x)
    a4 = sinpi(t - 2 * x)
    a5 = sinpi(2 * t - 4 * x)
    a6 = exp(t / 2)
    dq1 = -pi^2 * D^2 * (a4 + 2 * pi * a3) * a6 / 3 - 3 * pi * a2 * exp(t) * a5 / 2 +
          2 * pi * a1 * a6 * a3 - a6 * a4 / 2 - pi * a6 * a3

    return SVector(dq1)
end

function create_cache(mesh, equations::BBMEquation1D,
                      solver, initial_condition,
                      ::BoundaryConditionPeriodic,
                      RealT, uEltype)
    D = equations.D
    invImD2 = lu(I - 1 / 6 * D^2 * sparse(solver.D2))
    eta2 = zeros(RealT, nnodes(mesh))
    eta2_x = zero(eta2)
    eta_x = zero(eta2)
    etaeta_x = zero(eta2)
    eta_xx = zero(eta2)
    g = gravity_constant(equations)
    c_0 = sqrt(g * D)
    c_1 = sqrt(g / D)
    cache = (; invImD2, eta2, eta2_x, eta_x, etaeta_x, eta_xx, c_0, c_1,
             solver.D1, solver.D2)
    if solver.D1 isa PeriodicUpwindOperators
        eta_x_upwind = zero(eta2)
        cache = (; cache..., eta_x_upwind)
    end
    return cache
end

# Discretization that conserves the mass for eta and the modified energy for periodic boundary conditions, see
# - Hendrik Ranocha, Dimitrios Mitsotakis and David I. Ketcheson (2020)
#   A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations
#   [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)
function rhs!(dq, q, t, mesh, equations::BBMEquation1D, initial_condition,
              ::BoundaryConditionPeriodic, source_terms, solver, cache)
    (; invImD2, eta2, eta2_x, eta_x, etaeta_x, c_0, c_1) = cache
    if solver.D1 isa PeriodicUpwindOperators
        (; eta_x_upwind) = cache
    end

    eta, = q.x
    deta, = dq.x

    @trixi_timeit timer() "hyperbolic" begin
        @.. eta2 = eta^2
        if solver.D1 isa PeriodicDerivativeOperator ||
           solver.D1 isa UniformPeriodicCoupledOperator ||
           solver.D1 isa FourierDerivativeOperator
            mul!(eta2_x, solver.D1, eta2)
            mul!(eta_x, solver.D1, eta)
            if equations.split_form
                @.. etaeta_x = eta * eta_x
                @.. deta = -(0.5 * c_1 * (eta2_x + etaeta_x) + c_0 * eta_x)
            else
                @.. deta = -(0.75 * c_1 * eta2_x + c_0 * eta_x)
            end
        elseif solver.D1 isa PeriodicUpwindOperators
            mul!(eta2_x, solver.D1.central, eta2)
            mul!(eta_x_upwind, solver.D1.minus, eta)
            mul!(eta_x, solver.D1.central, eta)
            if equations.split_form
                @.. etaeta_x = eta * eta_x
                @.. deta = -(0.5 * c_1 * (eta2_x + etaeta_x) + c_0 * eta_x_upwind)
            else
                @.. deta = -(0.75 * c_1 * eta2_x + c_0 * eta_x_upwind)
            end
        else
            @error "unknown type of first-derivative operator: $(typeof(solver.D1))"
        end
    end

    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)

    @trixi_timeit timer() "elliptic" begin
        solve_system_matrix!(deta, invImD2, equations)
    end
    return nothing
end

"""
    energy_total_modified!(e, q_global, equations::BBMEquation1D, cache)

Return the modified total energy `e` of the primitive variables `q_global` for the
[`BBMEquation1D`](@ref). The `energy_total_modified`
is a conserved quantity (for periodic boundary conditions).

It is given by
```math
\\frac{1}{2} \\eta(\\eta - \\eta_{xx}).
```

`q_global` is a vector of the primitive variables at ALL nodes.
`cache` needs to hold the SBP operators used by the `solver`.

See also [`energy_total_modified`](@ref).
"""
function energy_total_modified!(e, q_global, equations::BBMEquation1D, cache)
    eta, = q_global.x

    (; D1, D2, eta_xx, tmp1) = cache
    if D1 isa PeriodicUpwindOperators
        mul!(tmp1, D1.minus, eta)
        mul!(eta_xx, D1.plus, tmp1)
    else
        mul!(eta_xx, D2, eta)
    end

    @.. e = 0.5 * eta * (eta - eta_xx)
    return e
end

"""
    hamiltonian!(H, q_global, equations::BBMEquation1D, cache)

Return the Hamiltonian `H` of the primitive variables `q_global` for the
[`BBMEquation1D`](@ref). The Hamiltonian is given by
```math
\\frac{1}{4}\\sqrt{\frac{g}{D}}\\eta^3 + \\frac{1}{2}\\sqrt{gD}\\eta^2.
```

`q_global` is a vector of the primitive variables at ALL nodes.

See also [`hamiltonian`](@ref).
"""
function hamiltonian!(H, q_global, equations::BBMEquation1D, cache)
    eta, = q_global.x
    (; c_0, c_1, eta2, tmp1) = cache
    @.. eta2 = c_0 * eta^2
    @.. tmp1 = c_1 * eta^3
    @.. H = 1 / 4 * tmp1 + 1 / 2 * eta2
    return H
end
