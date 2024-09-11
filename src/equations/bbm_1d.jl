# TODO: Use this nondimensionalized version or use dimensional version including gravity and depth?
@doc raw"""
    BBMEquation1D(; bathymetry_type = bathymetry_flat, eta0 = 0.0)

BBM (Benjamin–Bona–Mahony) equation in one spatial dimension. The equations are given by
```math
\begin{aligned}
  \eta_t + \eta\eta_x + \eta_x - \eta_{xxt} &= 0.
\end{aligned}
```

The unknown quantity of the BBM equation is the total water height ``\eta``. The water height
is given by ``h = \eta - \eta_0``, where ``\eta_0`` is the constant still-water surface.

Currently, the equations only support a flat bathymetry ``b = 0``, see [`bathymetry_flat`](@ref).

The BBM equation is first described in Benjamin, Bona, and Mahony (1972).
The semidiscretization implemented here is developed in Ranocha, Mitsotakis, and Ketcheson (2020).
It conserves
- the total water mass (integral of ``h``) as a linear invariant
- the total modified energy

for periodic boundary conditions.

- Thomas B. Benjamin, Jerry L. Bona and John J. Mahony (1972)
  Model equations for long waves in nonlinear dispersive systems
  [DOI: 10.1098/rsta.1972.0032](https://doi.org/10.1098/rsta.1972.0032)
- Hendrik Ranocha, Dimitrios Mitsotakis and David I. Ketcheson (2020)
  A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations
  [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)
"""
struct BBMEquation1D{Bathymetry <: AbstractBathymetry, RealT <: Real} <:
       AbstractBBMEquation{1, 1}
    bathymetry_type::Bathymetry # type of bathymetry
    eta0::RealT # constant still-water surface
end

function BBMEquation1D(; bathymetry_type = bathymetry_flat, eta0 = 0.0)
    BBMEquation1D(bathymetry_type, eta0)
end

"""
    initial_condition_convergence_test(x, t, equations::BBMEquation1D, mesh)

A travelling-wave solution used for convergence tests in a periodic domain.

See section 4.1.3 in (there is an error in paper: it should be sech^2 instead of cosh):
- Hendrik Ranocha, Dimitrios Mitsotakis and David I. Ketcheson (2020)
  A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations
  [DOI: 10.4208/cicp.OA-2020-0119](https://doi.org/10.4208/cicp.OA-2020-0119)
"""
function initial_condition_convergence_test(x, t, equations::BBMEquation1D, mesh)
    c = 1.2
    A = 3 * (c - 1)
    K = 0.5 * sqrt(1 - 1 / c)
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
    a3 = cospi(t - 2 * x)
    a4 = sinpi(t - 2 * x)
    a5 = sinpi(2 * t - 4 * x)
    a14 = exp(t / 2)
    dq1 = -2 * pi^2 * (a4 + 2 * pi * a3) * a14 - a14 * a4 / 2 + pi * a14 * a3 -
          pi * exp(t) * a5

    return SVector(dq1)
end

function create_cache(mesh, equations::BBMEquation1D,
                      solver, initial_condition,
                      ::BoundaryConditionPeriodic,
                      RealT, uEltype)
    invImD2 = lu(I - sparse(solver.D2))
    eta2 = zeros(RealT, nnodes(mesh))
    eta2_x = zero(eta2)
    eta_x = zero(eta2)
    etaeta_x = zero(eta2)
    cache = (; invImD2, eta2, eta2_x, eta_x, etaeta_x, solver.D1)
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
    (; invImD2, eta2, eta2_x, eta_x, etaeta_x) = cache
    if solver.D1 isa PeriodicUpwindOperators
        (; eta_x_upwind) = cache
    end

    eta, = q.x
    deta, = dq.x
    # energy and mass conservative semidiscretization
    @trixi_timeit timer() "hyperbolic" begin
        @.. eta2 = eta^2
        if solver.D1 isa PeriodicDerivativeOperator ||
           solver.D1 isa UniformPeriodicCoupledOperator ||
           solver.D1 isa FourierDerivativeOperator
            mul!(eta2_x, solver.D1, eta2)
            mul!(eta_x, solver.D1, eta)
            @.. etaeta_x = eta * eta_x
            # split form
            @.. deta = -(1 / 3 * (eta2_x + etaeta_x) + eta_x)
        elseif solver.D1 isa PeriodicUpwindOperators
            mul!(eta2_x, solver.D1.central, eta2)
            mul!(eta_x_upwind, solver.D1.minus, eta)
            mul!(eta_x, solver.D1.central, eta)
            @.. etaeta_x = eta * eta_x
            # split form
            @.. deta = -(1 / 3 * (eta2_x + etaeta_x) + eta_x_upwind)
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
    energy_total_modified!(e, q, equations::BBMEquation1D, cache)

Return the modified total energy `e` of the primitive variables `q_global` for the
[`BBMEquation1D`](@ref). The `energy_total_modified`
is a conserved quantity (for periodic boundary conditions).

It is given by
```math
\\frac{1}{2} \\eta^2 + \\frac{1}{2} \\eta_x^2.
```

`q_global` is a vector of the primitive variables at ALL nodes.
`cache` needs to hold the SBP operators used by the `solver`.

See also [`energy_total_modified`](@ref).
"""
function energy_total_modified!(e, q_global, equations::BBMEquation1D, cache)
    eta, = q_global.x

    (; D1, eta_x) = cache
    if D1 isa PeriodicUpwindOperators
        mul!(eta_x, D1.central, eta)
    else
        mul!(eta_x, D1, eta)
    end

    @.. e = 0.5 * (eta^2 + eta_x^2)
    return nothing
end

"""
    invariant_cubic(q, equations::BBMEquation1D)

Return the cubic invariant of the primitive variables `q` for the
[`BBMEquation1D`](@ref). The cubic invariant is given by
```math
(\\eta + 1)^3.
```
"""
function invariant_cubic(q, equations::BBMEquation1D)
    eta = waterheight_total(q, equations)
    return (eta + 1)^3
end

varnames(::typeof(invariant_cubic), equations) = ("(η + 1)³",)
