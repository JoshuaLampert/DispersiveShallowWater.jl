@doc raw"""
    SerreGreenNaghdiEquations1D(; gravity_constant, eta0 = 0.0)

Serre-Green-Naghdi system in one spatial dimension.
The equations for flat bathymetry are given by
```math
\begin{aligned}
  h_t + (h v)_x &= 0,\\
  h v_t - \frac{1}{3} (h^3 v_{tx})_x + \frac{1}{2} g (h^2)_x + \frac{1}{2} h (v^2)_x + p_x &= 0,\\
  p &= \frac{1}{3} h^3 v_{x}^2 - \frac{1}{3} h^3 v v_{xx}.
\end{aligned}
```
The unknown quantities of the Serre-Green-Naghdi equations are the
total water height ``\eta = h + b`` and the velocity ``v``.
The gravitational constant is denoted by `g` and the bottom topography
(bathymetry) ``b = \eta_0 - D``. The water height above the bathymetry
is therefore given by ``h = \eta - \eta_0 + D``. The Serre-Green-Naghdi
equations are only implemented for ``\eta_0 = 0``.
The total water height is therefore given by ``\eta = h + b``.

References for the Serre-Green-Naghdi system can be found in
- Serre (1953)
  Contribution â l'étude des écoulements permanents et variables dans les canaux
  [DOI: 10.1051/lhb/1953034](https://doi.org/10.1051/lhb/1953034)
- Green and Naghdi (1976)
  A derivation of equations for wave propagation in water of variable depth
  [DOI: 10.1017/S0022112076002425](https://doi.org/10.1017/S0022112076002425)

The semidiscretization implemented here conserves
- the total water mass (integral of h) as a linear invariant
- the total momentum (integral of h v) as a nonlinear invariant
- the total energy

for periodic boundary conditions, see
- Hendrik Ranocha and Mario Ricchiuto (2024)
  Structure-preserving approximations of the Serre-Green-Naghdi
  equations in standard and hyperbolic form
  [arXiv: 2408.02665](https://arxiv.org/abs/2408.02665)
"""
struct SerreGreenNaghdiEquations1D{Bathymetry <: AbstractBathymetry, RealT <: Real} <:
       AbstractSerreGreenNaghdiEquations{1, 3}
    bathymetry::Bathymetry # type of bathymetry
    gravity::RealT # gravitational constant
    eta0::RealT # constant still-water surface
end

function SerreGreenNaghdiEquations1D(; gravity_constant, eta0 = 0.0)
    eta0 == 0.0 ||
        @warn "The still-water surface needs to be 0 for the Serre-Green-Naghdi equations"
    SerreGreenNaghdiEquations1D(bathymetry_flat, gravity_constant, eta0)
end

varnames(::typeof(prim2prim), ::SerreGreenNaghdiEquations1D) = ("η", "v", "D")
varnames(::typeof(prim2cons), ::SerreGreenNaghdiEquations1D) = ("h", "hv", "b")

"""
    initial_condition_convergence_test(x, t, equations::SerreGreenNaghdiEquations1D, mesh)

A soliton solution used for convergence tests in a periodic domain.
"""
function initial_condition_convergence_test(x, t, equations::SerreGreenNaghdiEquations1D,
                                            mesh)
    g = equations.gravity

    # setup parameters data
    h1 = 1.0
    h2 = 1.2
    c = sqrt(g * h2)

    x_t = mod(x - c * t - xmin(mesh), xmax(mesh) - xmin(mesh)) + xmin(mesh)

    h = h1 + (h2 - h1) * sech(x_t / 2 * sqrt(3 * (h2 - h1) / (h1^2 * h2)))^2
    v = c * (1 - h1 / h)

    return SVector(h, v, 0)
end

function create_cache(mesh,
                      equations::SerreGreenNaghdiEquations1D{BathymetryFlat},
                      solver,
                      initial_condition,
                      ::BoundaryConditionPeriodic,
                      RealT, uEltype)
    D = solver.D1

    # create temporary storage
    h = ones(RealT, nnodes(mesh))
    h_x = zero(h)
    v_x = zero(h)
    h2_x = zero(h)
    hv_x = zero(h)
    v2_x = zero(h)
    h2_v_vx_x = zero(h)
    h_vx_x = zero(h)
    p_x = zero(h)
    tmp = zero(h)
    M_h = zero(h)
    M_h3_3 = zero(h)

    if D isa PeriodicUpwindOperators
        v_x_upwind = zero(h)

        Dmat_minus = sparse(D.minus)

        # Floating point errors accumulate a bit and the system matrix
        # is not necessarily perfectly symmetric but only up to
        # round-off errors. We wrap it here to avoid issues with the
        # factorization.
        @. M_h = h
        scale_by_mass_matrix!(M_h, D)
        @. M_h3_3 = (1 / 3) * h^3
        scale_by_mass_matrix!(M_h3_3, D)
        system_matrix = Symmetric(Diagonal(M_h)
                                  +
                                  Dmat_minus' * Diagonal(M_h3_3) * Dmat_minus)
        factorization = cholesky(system_matrix)

        cache = (; h_x, v_x, v_x_upwind, h2_x, hv_x, v2_x,
                 h2_v_vx_x, h_vx_x, p_x, tmp,
                 M_h, M_h3_3,
                 D, Dmat_minus, factorization)
    else
        if D isa FourierDerivativeOperator
            Dmat = Matrix(D)

            cache = (; h_x, v_x, h2_x, hv_x, v2_x,
                     h2_v_vx_x, h_vx_x, p_x, tmp,
                     M_h, M_h3_3,
                     D, Dmat)
        else
            Dmat = sparse(D)

            # Floating point errors accumulate a bit and the system matrix
            # is not necessarily perfectly symmetric but only up to
            # round-off errors. We wrap it here to avoid issues with the
            # factorization.
            @. M_h = h
            scale_by_mass_matrix!(M_h, D)
            @. M_h3_3 = (1 / 3) * h^3
            scale_by_mass_matrix!(M_h3_3, D)
            system_matrix = Symmetric(Diagonal(M_h)
                                      +
                                      Dmat' * Diagonal(M_h3_3) * Dmat)
            factorization = cholesky(system_matrix)

            cache = (; h_x, v_x, h2_x, hv_x, v2_x,
                     h2_v_vx_x, h_vx_x, p_x, tmp,
                     M_h, M_h3_3,
                     D, Dmat, factorization)
        end
    end

    return cache
end

# Discretization that conserves
# - the total water mass (integral of h) as a linear invariant
# - the total momentum (integral of h v) as a nonlinear invariant
# - the total energy
# for periodic boundary conditions, see
# - Hendrik Ranocha and Mario Ricchiuto (2024)
#   Structure-preserving approximations of the Serre-Green-Naghdi
#   equations in standard and hyperbolic form
#   [arXiv: 2408.02665](https://arxiv.org/abs/2408.02665)
# TODO: Implement source terms
# TODO: Implement variable bathymetry
function rhs!(dq, q, t, mesh,
              equations::SerreGreenNaghdiEquations1D{BathymetryFlat},
              initial_condition,
              ::BoundaryConditionPeriodic,
              source_terms::Nothing,
              solver, cache)
    if cache.D isa PeriodicUpwindOperators
        rhs_sgn_flat_upwind!(dq, q, equations, source_terms, cache)
    else
        rhs_sgn_flat_central!(dq, q, equations, source_terms, cache)
    end

    return nothing
end

function rhs_sgn_flat_central!(dq, q, equations, source_terms, cache)
    # Unpack physical parameters and SBP operator `D` as well as the
    # SBP operator in sparse matrix form `Dmat`
    g = equations.gravity
    (; D, Dmat) = cache

    # `q` and `dq` are `ArrayPartition`s. They collect the individual
    # arrays for the water height `h` and the velocity `v`.
    h, v = q.x
    dh, dv = dq.x

    @trixi_timeit timer() "hyperbolic terms" begin
        # Compute all derivatives required below
        (; h_x, v_x, h2_x, hv_x, v2_x, h2_v_vx_x,
        h_vx_x, p_x, tmp, M_h, M_h3_3) = cache

        mul!(h_x, D, h)
        mul!(v_x, D, v)
        @. tmp = h^2
        mul!(h2_x, D, tmp)
        @. tmp = h * v
        mul!(hv_x, D, tmp)
        @. tmp = v^2
        mul!(v2_x, D, tmp)

        @. tmp = h^2 * v * v_x
        mul!(h2_v_vx_x, D, tmp)
        @. tmp = h * v_x
        mul!(h_vx_x, D, tmp)
        inv6 = 1 / 6
        @. tmp = (0.5 * h^2 * (h * v_x + h_x * v) * v_x
                  -
                  inv6 * h * h2_v_vx_x
                  -
                  inv6 * h^2 * v * h_vx_x)
        mul!(p_x, D, tmp)

        # Plain: h_t + (h v)_x = 0
        #
        # Split form for energy conservation:
        # h_t + h_x v + h v_x = 0
        @. dh = -(h_x * v + h * v_x)

        # Plain: h v_t + ... = 0
        #
        # Split form for energy conservation:
        @. tmp = -(g * h2_x - g * h * h_x
                   +
                   0.5 * h * v2_x
                   -
                   0.5 * v^2 * h_x
                   +
                   0.5 * hv_x * v
                   -
                   0.5 * h * v * v_x
                   +
                   p_x)
    end

    @trixi_timeit timer() "assembling elliptic operator" begin
        # The code below is equivalent to
        #   dv .= (Diagonal(h) - Dmat * Diagonal(1/3 .* h.^3) * Dmat) \ tmp
        # but faster since the symbolic factorization is reused.
        # Floating point errors accumulate a bit and the system matrix
        # is not necessarily perfectly symmetric but only up to round-off
        # errors. We wrap it here to avoid issues with the factorization.
        @. M_h = h
        scale_by_mass_matrix!(M_h, D)
        inv3 = 1 / 3
        @. M_h3_3 = inv3 * h^3
        scale_by_mass_matrix!(M_h3_3, D)
        system_matrix = Symmetric(Diagonal(M_h)
                                  +
                                  Dmat' * Diagonal(M_h3_3) * Dmat)
    end

    @trixi_timeit timer() "solving elliptic system" begin
        if issparse(system_matrix)
            (; factorization) = cache
            cholesky!(factorization, system_matrix; check = false)
            if issuccess(factorization)
                scale_by_mass_matrix!(tmp, D)
                dv .= factorization \ tmp
            else
                # The factorization may fail if the time step is too large
                # and h becomes negative.
                fill!(dv, NaN)
            end
        else
            factorization = cholesky!(system_matrix)
            scale_by_mass_matrix!(tmp, D)
            ldiv!(dv, factorization, tmp)
        end
    end

    return nothing
end

function rhs_sgn_flat_upwind!(dq, q, equations, source_terms, cache)
    # Unpack physical parameters and SBP operator `D` as well as the
    # SBP upwind operator in sparse matrix form `Dmat_minus`
    g = equations.gravity
    (; Dmat_minus) = cache
    D_upwind = cache.D
    D = D_upwind.central

    # `q` and `dq` are `ArrayPartition`s. They collect the individual
    # arrays for the water height `h` and the velocity `v`.
    h, v = q.x
    dh, dv = dq.x

    @trixi_timeit timer() "hyperbolic terms" begin
        # Compute all derivatives required below
        (; h_x, v_x, v_x_upwind, h2_x, hv_x, v2_x,
        h2_v_vx_x, h_vx_x, p_x, tmp,
        M_h, M_h3_3) = cache

        mul!(h_x, D, h)
        mul!(v_x, D, v)
        mul!(v_x_upwind, D_upwind.minus, v)
        @. tmp = h^2
        mul!(h2_x, D, tmp)
        @. tmp = h * v
        mul!(hv_x, D, tmp)
        @. tmp = v^2
        mul!(v2_x, D, tmp)

        @. tmp = h^2 * v * v_x
        mul!(h2_v_vx_x, D, tmp)
        @. tmp = h * v_x
        mul!(h_vx_x, D, tmp)
        # p_+
        @. tmp = 0.5 * h^2 * (h * v_x + h_x * v) * v_x_upwind
        mul!(p_x, D_upwind.plus, tmp)
        # p_0
        minv6 = -1 / 6
        @. tmp = minv6 * (h * h2_v_vx_x
                          +
                          h^2 * v * h_vx_x)
        mul!(p_x, D, tmp, 1.0, 1.0)

        # Plain: h_t + (h v)_x = 0
        #
        # Split form for energy conservation:
        # h_t + h_x v + h v_x = 0
        @. dh = -(h_x * v + h * v_x)

        # Plain: h v_t + ... = 0
        #
        # Split form for energy conservation:
        @. tmp = -(g * h2_x - g * h * h_x
                   +
                   0.5 * h * v2_x
                   -
                   0.5 * v^2 * h_x
                   +
                   0.5 * hv_x * v
                   -
                   0.5 * h * v * v_x
                   +
                   p_x)
    end

    @trixi_timeit timer() "assembling elliptic operator" begin
        # The code below is equivalent to
        #   dv .= (Diagonal(h) - Dmat_plus * Diagonal(1/3 .* h.^3) * Dmat_minus) \ tmp
        # but faster since the symbolic factorization is reused.
        # Floating point errors accumulate a bit and the system matrix
        # is not necessarily perfectly symmetric but only up to round-off errors.
        # We wrap it here to avoid issues with the factorization.
        @. M_h = h
        scale_by_mass_matrix!(M_h, D)
        inv3 = 1 / 3
        @. M_h3_3 = inv3 * h^3
        scale_by_mass_matrix!(M_h3_3, D)
        system_matrix = Symmetric(Diagonal(M_h)
                                  +
                                  Dmat_minus' * Diagonal(M_h3_3) * Dmat_minus)
    end

    @trixi_timeit timer() "solving elliptic system" begin
        (; factorization) = cache
        cholesky!(factorization, system_matrix; check = false)
        if issuccess(factorization)
            scale_by_mass_matrix!(tmp, D)
            dv .= factorization \ tmp
        else
            # The factorization may fail if the time step is too large
            # and h becomes negative.
            fill!(dv, NaN)
        end
    end

    return nothing
end

@inline function prim2cons(q, equations::SerreGreenNaghdiEquations1D)
    h = waterheight(q, equations)
    v = velocity(q, equations)
    b = bathymetry(q, equations)

    hv = h * v
    return SVector(h, hv, b)
end

@inline function cons2prim(u, equations::SerreGreenNaghdiEquations1D)
    h, hv, b = u

    eta = h + b
    v = hv / h
    D = equations.eta0 - b
    return SVector(eta, v, D)
end

@inline function waterheight_total(q, equations::SerreGreenNaghdiEquations1D)
    return q[1]
end

@inline function velocity(q, equations::SerreGreenNaghdiEquations1D)
    return q[2]
end

@inline function bathymetry(q, equations::SerreGreenNaghdiEquations1D)
    D = q[3]
    return equations.eta0 - D
end

@inline function waterheight(q, equations::SerreGreenNaghdiEquations1D)
    return waterheight_total(q, equations) - bathymetry(q, equations)
end

# The entropy/energy takes the whole `q` for every point in space
"""
    energy_total_modified(q_global, equations::SerreGreenNaghdiEquations1D, cache)

Return the modified total energy of the primitive variables `q_global` for the
[`SerreGreenNaghdiEquations1D`](@ref).
It contains an additional term containing a
derivative compared to the usual [`energy_total`](@ref) modeling
non-hydrostatic contributions. The [`energy_total_modified`](@ref)
is a conserved quantity (for periodic boundary conditions).

For a [`bathymetry_flat`](@ref) the total energy is given by
```math
\\frac{1}{2} g h^2 + \\frac{1}{2} h v^2 + \\frac{1}{6} h^3 v_x^2.
```

`q_global` is a vector of the primitive variables at ALL nodes.
`cache` needs to hold the SBP operators used by the `solver`.
"""
function energy_total_modified(q_global,
                               equations::SerreGreenNaghdiEquations1D,
                               cache)
    # unpack physical parameters and SBP operator `D`
    g = equations.gravity
    (; D, v_x) = cache

    # `q_global` is an `ArrayPartition`. It collects the individual arrays for
    # the water height `h` and the velocity `v`.
    h, v = q_global.x

    N = length(v)
    e = zeros(eltype(q), N)

    # 1/2 g h^2 + 1/2 h v^2 + 1/6 h^3 v_x^2
    if D isa PeriodicUpwindOperators
        mul!(v_x, D.minus, v)
    else
        mul!(v_x, D, v)
    end

    @. e = 1 / 2 * g * h^2 + 1 / 2 * h * v^2 + 1 / 6 * h^3 * v_x^2

    return e
end
