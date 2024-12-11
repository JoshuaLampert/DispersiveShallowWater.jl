"""
    AbstractEquations{NDIMS, NVARS}

An abstract supertype of specific equations such as the BBM-BBM equations.
The type parameters encode the number of spatial dimensions (`NDIMS`) and the
number of primary variables (`NVARS`) of the physics model.

See also [`AbstractShallowWaterEquations`](@ref).
"""
abstract type AbstractEquations{NDIMS, NVARS} end

# Retrieve number of variables from equation instance
@inline nvariables(::AbstractEquations{NDIMS, NVARS}) where {NDIMS, NVARS} = NVARS

"""
    eachvariable(equations::AbstractEquations)

Return an iterator over the indices that specify the location in relevant data structures
for the variables in `equations`. In particular, not the variables themselves are returned.
"""
@inline eachvariable(equations::AbstractEquations) = Base.OneTo(nvariables(equations))

@doc raw"""
    AbstractShallowWaterEquations{NDIMS, NVARS}

An abstract supertype of all equation system that contain the classical
shallow water equations as a subsystem, e.g., the
[`BBMBBMEquations1D`](@ref), the [`SvaerdKalischEquations1D`](@ref),
and the [`SerreGreenNaghdiEquations1D`](@ref).
In 1D, the shallow water equations with flat bathymetry are given by
```math
\begin{aligned}
  h_t + (h v)_x &= 0,\\
  h v_t + \frac{1}{2} g (h^2)_x + \frac{1}{2} h (v^2)_x &= 0,
\end{aligned}
```
where ``h`` is the [`waterheight`](@ref),
``v`` the [`velocity`](@ref), and
``g`` the [`gravity_constant`](@ref).
"""
abstract type AbstractShallowWaterEquations{NDIMS, NVARS} <: AbstractEquations{NDIMS, NVARS} end

"""
    get_name(equations::AbstractEquations)

Return the canonical, human-readable name for the given system of equations.
# Examples
```jldoctest
julia> DispersiveShallowWater.get_name(BBMBBMEquations1D(gravity_constant=1.0))
"BBMBBMEquations1D-BathymetryVariable"
```
"""
get_name(equations::AbstractEquations) = equations |> typeof |> nameof |> string

function get_name(equations::AbstractShallowWaterEquations)
    name = equations |> typeof |> nameof |> string
    bathymetry_type = equations.bathymetry_type
    if bathymetry_type isa BathymetryFlat
        return name
    else # variable bathymetry_type
        return name * "-" * string(nameof(typeof(bathymetry_type)))
    end
end

"""
    varnames(conversion_function, equations)

Return the list of variable names when applying `conversion_function` to the
primitive variables associated to `equations`.
Common choices of the `conversion_function` are [`prim2prim`](@ref),
[`prim2cons`](@ref), and [`prim2phys`](@ref).
"""
function varnames end

"""
    prim2prim(q, equations)

Return the primitive variables `q`. While this function is as trivial as `identity`,
it is also as useful.
"""
@inline prim2prim(q, ::AbstractEquations) = q

varnames(::typeof(prim2prim), ::AbstractShallowWaterEquations) = ("η", "v", "D")
varnames(::typeof(prim2prim), ::AbstractEquations{1, 1}) = ("η",)

"""
    prim2cons(q, equations)

Convert the primitive variables `q` to the conserved variables for a given set of
`equations`. `q` is a vector type of the correct length `nvariables(equations)`.
Notice the function doesn't include any error checks for the purpose of efficiency,
so please make sure your input is correct.
The inverse conversion is performed by [`cons2prim`](@ref).
"""
@inline function prim2cons(q, equations::AbstractShallowWaterEquations)
    h = waterheight(q, equations)
    v = velocity(q, equations)
    b = bathymetry(q, equations)

    hv = h * v
    return SVector(h, hv, b)
end

function prim2cons(q, equations::AbstractEquations{1, 1})
    h = waterheight(q, equations)
    return SVector(h)
end

varnames(::typeof(prim2cons), ::AbstractShallowWaterEquations) = ("h", "hv", "b")
varnames(::typeof(prim2cons), ::AbstractEquations{1, 1}) = ("h",)

"""
    cons2prim(u, equations)

Convert the conserved variables `u` to the primitive variables for a given set of
`equations`. `u` is a vector type of the correct length `nvariables(equations)`.
Notice the function doesn't include any error checks for the purpose of efficiency,
so please make sure your input is correct.
The inverse conversion is performed by [`prim2cons`](@ref).
"""
@inline function cons2prim(u, equations::AbstractShallowWaterEquations)
    h, hv, b = u

    eta = h + b
    v = hv / h
    eta0 = equations.eta0
    D = eta0 - b
    return SVector(eta, v, D)
end

function cons2prim(u, equations::AbstractEquations{1, 1})
    h, = u
    eta0 = equations.eta0
    D = equations.D
    b = eta0 - D
    eta = h + b
    return SVector(eta)
end

"""
    waterheight_total(q, equations)

Return the total waterheight of the primitive variables `q` for a given set of
`equations`, i.e., the [`waterheight`](@ref) ``h`` plus the
[`bathymetry`](@ref) ``b``.

`q` is a vector of the primitive variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.
"""
@inline function waterheight_total(q, ::AbstractEquations)
    return q[1]
end

varnames(::typeof(waterheight_total), equations) = ("η",)

"""
    waterheight(q, equations)

Return the waterheight of the primitive variables `q` for a given set of
`equations`, i.e., the waterheight ``h`` above the bathymetry ``b``.

`q` is a vector of the primitive variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.

See also [`waterheight_total`](@ref), [`bathymetry`](@ref).
"""
@inline function waterheight(q, equations::AbstractEquations)
    return waterheight_total(q, equations) - bathymetry(q, equations)
end

varnames(::typeof(waterheight), equations) = ("h",)

"""
    velocity(q, equations)

Return the velocity of the primitive variables `q` for a given set of
`equations`.

`q` is a vector of the primitive variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.
"""
@inline function velocity(q, ::AbstractShallowWaterEquations)
    return q[2]
end

varnames(::typeof(velocity), equations) = ("v",)

"""
    momentum(q, equations)

Return the momentum/discharge of the primitive variables `q` for a given set of
`equations`, i.e., the [`waterheight`](@ref) times the [`velocity`](@ref).

`q` is a vector of the primitive variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.
"""
@inline function momentum(q, equations::AbstractEquations)
    return waterheight(q, equations) * velocity(q, equations)
end

varnames(::typeof(momentum), equations) = ("P",)

"""
    discharge(q, equations)

See [`momentum`](@ref).
"""
@inline discharge(q, equations::AbstractEquations) = momentum(q, equations)

varnames(::typeof(discharge), equations) = ("P",)

"""
    still_water_surface(q, equations)

Return the still water surface ``\\eta_0`` (lake at rest)
for a given set of `equations`.
"""
@inline function still_water_surface(q, equations::AbstractEquations)
    return equations.eta0
end

@inline function still_waterdepth(q, ::AbstractShallowWaterEquations)
    return q[3]
end

still_waterdepth(q, equations::AbstractEquations{1, 1}) = equations.D

"""
    bathymetry(q, equations)

Return the bathymetry of the primitive variables `q` for a given set of
`equations`.

`q` is a vector of the primitive variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.
"""
@inline function bathymetry(q, equations::AbstractEquations)
    D = still_waterdepth(q, equations)
    eta0 = still_water_surface(q, equations)
    return eta0 - D
end

"""
    lake_at_rest_error(q, equations)

Calculate the error for the "lake-at-rest" test case where the
[`waterheight_total`](@ref) ``\\eta = h + b`` should
be a constant value over time (given by the value ``\\eta_0`` passed to the
`equations` when constructing them).
"""
@inline function lake_at_rest_error(q, equations::AbstractEquations)
    eta = waterheight_total(q, equations)
    eta0 = still_water_surface(q, equations)
    return abs(eta0 - eta)
end

"""
    gravity_constant(equations)

Return the gravity constant ``g`` for a given set of `equations`.
"""
@inline function gravity_constant(equations::AbstractEquations)
    return equations.gravity
end

"""
    energy_total(q, equations)

Return the total energy of the primitive variables `q` for a given set of
`equations`. For all [`AbstractShallowWaterEquations`](@ref), the total
energy is given by the sum of the kinetic and potential energy of the
shallow water subsystem, i.e.,
```math
\\frac{1}{2} h v^2 + \\frac{1}{2} g \\eta^2
```
in 1D, where ``h`` is the [`waterheight`](@ref),
``\\eta = h + b`` the [`waterheight_total`](@ref),
``v`` the [`velocity`](@ref), and ``g`` the [`gravity_constant`](@ref).

`q` is a vector of the primitive variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.
"""
@inline function energy_total(q, equations::AbstractShallowWaterEquations)
    h = waterheight(q, equations)
    eta = waterheight_total(q, equations)
    v = velocity(q, equations)
    return 0.5f0 * h * v^2 + 0.5f0 * gravity_constant(equations) * eta^2
end

varnames(::typeof(energy_total), equations) = ("e_total",)

"""
    entropy(q, equations)

Return the entropy of the primitive variables `q` for a given set of
`equations`. For all [`AbstractShallowWaterEquations`](@ref), the `entropy`
is just the [`energy_total`](@ref).

`q` is a vector of the primitive variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.
"""
function entropy(q, equations::AbstractShallowWaterEquations)
    return energy_total(q, equations)
end

varnames(::typeof(entropy), equations) = ("U",)

# The modified entropy/total energy takes the whole `q_global` for every point in space
"""
    energy_total_modified(q_global, equations, cache)

Return the modified total energy of the primitive variables `q_global` for the
`equations`. This modified total energy is a conserved quantity and can
contain additional terms compared to the usual [`energy_total`](@ref).
For example, for the [`SvaerdKalischEquations1D`](@ref) and the
[`SerreGreenNaghdiEquations1D`](@ref), it contains additional terms
depending on the derivative of the velocity ``v_x`` modeling non-hydrostatic
contributions. For equations which are not `AbstractShallowWaterEquations`,
the modified total energy does not have to be an extension of the usual
[`energy_total`](@ref) and does not have to be related to a physical energy.
However, it is still a conserved quantity for appropriate boundary conditions
and may contain derivatives of the solution.

`q_global` is a vector of the primitive variables at ALL nodes.
`cache` needs to hold the SBP operators used by the `solver` if non-hydrostatic
terms are present.

Internally, this function allocates a vector for the output and
calls [`DispersiveShallowWater.energy_total_modified!`](@ref).
"""
function energy_total_modified(q_global, equations::AbstractEquations, cache)
    e = similar(q_global.x[begin])
    return energy_total_modified!(e, q_global, equations, cache)
end

"""
    energy_total_modified!(e, q_global, equations, cache)

In-place version of [`energy_total_modified`](@ref).
"""
function energy_total_modified!(e, q_global, equations::AbstractEquations,
                                cache)
    # `q_global` is an `ArrayPartition` of the primitive variables at all nodes
    @assert nvariables(equations) == length(q_global.x)

    for i in eachindex(q_global.x[begin])
        e[i] = energy_total(get_node_vars(q_global, equations, i), equations)
    end

    return e
end

varnames(::typeof(energy_total_modified), equations) = ("e_modified",)

"""
    entropy_modified(q_global, equations, cache)

Alias for [`energy_total_modified`](@ref).
"""
@inline function entropy_modified(q_global, equations::AbstractEquations, cache)
    e = similar(q_global.x[begin])
    return entropy_modified!(e, q_global, equations, cache)
end

"""
    entropy_modified!(e, q_global, equations, cache)

In-place version of [`entropy_modified`](@ref).
"""
@inline entropy_modified!(e, q_global, equations, cache) = energy_total_modified!(e,
                                                                                  q_global,
                                                                                  equations,
                                                                                  cache)

varnames(::typeof(entropy_modified), equations) = ("U_modified",)

# The Hamiltonian takes the whole `q_global` for every point in space
"""
    hamiltonian(q_global, equations, cache)

Return the Hamiltonian of the primitive variables `q_global` for the
`equations`. The Hamiltonian is a conserved quantity and may contain
derivatives of the solution.

`q_global` is a vector of the primitive variables at ALL nodes.
`cache` needs to hold the SBP operators used by the `solver` if non-hydrostatic
terms are present.

Internally, this function allocates a vector for the output and
calls [`DispersiveShallowWater.hamiltonian!`](@ref).

!!! note
    This function is not necessarily implemented for all equations.
    See [`DispersiveShallowWater.hamiltonian!`](@ref) for further details
    of equations supporting it.
"""
function hamiltonian(q_global, equations::AbstractEquations, cache)
    H = similar(q_global.x[begin])
    return hamiltonian!(H, q_global, equations, cache)
end

varnames(::typeof(hamiltonian), equations) = ("H",)

# Add methods to show some information on systems of equations.
function Base.show(io::IO, equations::AbstractEquations)
    # Since this is not performance-critical, we can use `@nospecialize` to reduce latency.
    @nospecialize equations # reduce precompilation time

    print(io, get_name(equations), " with ")
    if nvariables(equations) == 1
        println(io, "one variable")
    else
        println(io, nvariables(equations), " variables")
    end
end

function Base.show(io::IO, ::MIME"text/plain", equations::AbstractEquations)
    # Since this is not performance-critical, we can use `@nospecialize` to reduce latency.
    @nospecialize equations # reduce precompilation time

    if get(io, :compact, false)
        show(io, equations)
    else
        println(io, get_name(equations))
        println(io, "#variables: ", nvariables(equations))
        for variable in eachvariable(equations)
            println("    variable " * string(variable), ": ",
                    varnames(prim2prim, equations)[variable])
        end
    end
end

@inline Base.ndims(::AbstractEquations{NDIMS}) where {NDIMS} = NDIMS

"""
    default_analysis_errors(equations)

Default analysis errors used by the [`AnalysisCallback`](@ref).
"""
default_analysis_errors(::AbstractEquations) = (:l2_error, :linf_error)

"""
    default_analysis_integrals(equations)

Default analysis integrals used by the [`AnalysisCallback`](@ref).
"""
default_analysis_integrals(::AbstractEquations) = Symbol[]

"""
    DispersiveShallowWater.is_hyperbolic_appproximation(equations)

Returns `Val{true}()` if the equations are a hyperbolic approximation
of another set of equations and `Val{false}()` otherwise (default).
For example, the [`HyperbolicSerreGreenNaghdiEquations1D`](@ref) are
a hyperbolic approximation of the [`SerreGreenNaghdiEquations1D`](@ref).

See also [`hyperbolic_approximation_limit`](@ref) and [`prim2phys`](@ref).

!!! note "Implementation details"
    This function is mostly used for some internal dispatch. For example,
    it allows to return a reduced set of variables from initial conditions
    for hyperbolic approximations.
"""
is_hyperbolic_appproximation(::AbstractEquations) = Val{false}()

"""
    DispersiveShallowWater.hyperbolic_approximation_limit(equations)

If the equations are a hyperbolic approximation of another set of equations,
return the equations of the limit system. Otherwise, return the input equations.

See also [`is_hyperbolic_appproximation`](@ref) and [`prim2phys`](@ref).

!!! note "Implementation details"
    This function is mostly used for some internal dispatch. For example,
    it allows to return a reduced set of variables from initial conditions
    for hyperbolic approximations.
"""
hyperbolic_approximation_limit(equations::AbstractEquations) = equations

"""
    prim2phys(q, equations)

Convert the primitive variables `q` to the physically meaningful variables
for a given set of `equations`. By default, this is the same as
[`prim2prim`](@ref) for most equations. However, some equations like the
[`HyperbolicSerreGreenNaghdiEquations1D`](@ref) return a reduced set of
variables since they are a hyperbolic approximation of another set of
equations (in this case the [`SerreGreenNaghdiEquations1D`](@ref)).

See also [`is_hyperbolic_appproximation`](@ref) and
[`hyperbolic_approximation_limit`](@ref).

`q` is a vector type of the correct length `nvariables(equations)`.
Notice the function doesn't include any error checks for the purpose of
efficiency, so please make sure your input is correct.
"""
prim2phys(q, equations::AbstractEquations) = prim2prim(q, equations)

function varnames(::typeof(prim2phys), equations::AbstractEquations)
    return varnames(prim2prim, equations)
end

abstract type AbstractBathymetry end

struct BathymetryFlat <: AbstractBathymetry end
"""
    bathymetry_flat = DispersiveShallowWater.BathymetryFlat()

A singleton struct indicating a flat bathymetry.

See also [`bathymetry_mild_slope`](@ref) and [`bathymetry_variable`](@ref).
"""
const bathymetry_flat = BathymetryFlat()

struct BathymetryMildSlope <: AbstractBathymetry end
"""
    bathymetry_mild_slope = DispersiveShallowWater.BathymetryMildSlope()

A singleton struct indicating a variable bathymetry with mild-slope
approximation. Typically, this means that some terms like ``b_x^2``
are neglected.

See also [`bathymetry_flat`](@ref) and [`bathymetry_variable`](@ref).
"""
const bathymetry_mild_slope = BathymetryMildSlope()

struct BathymetryVariable <: AbstractBathymetry end
"""
    bathymetry_variable = DispersiveShallowWater.BathymetryVariable()

A singleton struct indicating a variable bathymetry (without
mild-slope approximation).

See also [`bathymetry_flat`](@ref) and [`bathymetry_mild_slope`](@ref).
"""
const bathymetry_variable = BathymetryVariable()

# BBM equation
abstract type AbstractBBMEquation{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end
include("bbm_1d.jl")

# BBM-BBM equations
abstract type AbstractBBMBBMEquations{NDIMS, NVARS} <:
              AbstractShallowWaterEquations{NDIMS, NVARS} end
include("bbm_bbm_1d.jl")

# Svärd-Kalisch equations
abstract type AbstractSvaerdKalischEquations{NDIMS, NVARS} <:
              AbstractShallowWaterEquations{NDIMS, NVARS} end
include("svaerd_kalisch_1d.jl")

# Serre-Green-Naghdi equations
abstract type AbstractSerreGreenNaghdiEquations{NDIMS, NVARS} <:
              AbstractShallowWaterEquations{NDIMS, NVARS} end
include("serre_green_naghdi_1d.jl")
include("hyperbolic_serre_green_naghdi_1d.jl")

function solve_system_matrix!(dv, system_matrix, rhs,
                              ::Union{SvaerdKalischEquations1D,
                                      SerreGreenNaghdiEquations1D},
                              D1, cache)
    if issparse(system_matrix)
        (; factorization) = cache
        cholesky!(factorization, system_matrix; check = false)
        if issuccess(factorization)
            scale_by_mass_matrix!(rhs, D1)
            # see https://github.com/JoshuaLampert/DispersiveShallowWater.jl/issues/122
            dv .= system_matrix \ rhs[2:(end - 1)]
        else
            # The factorization may fail if the time step is too large
            # and h becomes negative.
            fill!(dv, NaN)
        end
    else
        factorization = cholesky!(system_matrix)
        scale_by_mass_matrix!(rhs, D1)
        ldiv!(dv, factorization, rhs)
    end
    return nothing
end

# Svärd-Kalisch equations for reflecting boundary conditions uses LU factorization
function solve_system_matrix_lu!(dv, system_matrix_lu, ::SvaerdKalischEquations1D, cache)
    (; factorization_lu) = cache
    lu!(factorization_lu, system_matrix_lu)
    ldiv!(factorization_lu, dv)
end

function solve_system_matrix!(dv, system_matrix, ::Union{BBMEquation1D, BBMBBMEquations1D})
    ldiv!(system_matrix, dv)
end

"""
    initial_condition_dingemans(x, t, equations::AbstractShallowWaterEquations, mesh)

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
function initial_condition_dingemans(x, t, equations::AbstractShallowWaterEquations, mesh)
    g = gravity_constant(equations)
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
    eta0 = equations.eta0
    eta = h + h0
    D = eta0 - b
    return SVector(eta, v, D)
end

"""
    initial_condition_discontinuous_well_balancedness(x, t, equations::AbstractShallowWaterEquations, mesh)

Setup a truly discontinuous bottom topography function for this academic
lake-at-rest test case of well-balancedness, i.e. `eta` is constant and `v` is
zero everywhere. The error for this lake-at-rest test case
`∫|η-η₀|` should be around machine roundoff.
"""
function initial_condition_discontinuous_well_balancedness(x, t,
                                                           equations::AbstractShallowWaterEquations,
                                                           mesh)
    eta0 = equations.eta0
    # Set the background values
    eta = eta0
    v = 0.0
    D = eta0 - 1.0

    # Setup a discontinuous bottom topography
    if x >= 0.5 && x <= 0.75
        D = eta0 - 1.5 - 0.5 * sinpi(2.0 * x)
    end

    return SVector(eta, v, D)
end
