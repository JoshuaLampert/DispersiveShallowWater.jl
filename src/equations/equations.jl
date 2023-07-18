"""
    AbstractEquations{NDIMS, NVARS}

An abstract supertype of specific equations such as the BBM-BBM equations.
The type parameters encode the number of spatial dimensions (`NDIMS`) and the
number of primary variables (`NVARS`) of the physics model.
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

"""
    get_name(equations::AbstractEquations)

Returns the canonical, human-readable name for the given system of equations.
# Examples
```jldoctest
julia> DispersiveShallowWater.get_name(BBMBBMEquations1D(gravity_constant=1.0))
"BBMBBMEquations1D"
```
"""
get_name(equations::AbstractEquations) = equations |> typeof |> nameof |> string

"""
    varnames(equations)

Return the list of variable names of `equations`.
"""
function varnames end

"""
    waterheight_total(u, equations)

Return the total waterheight of the conserved variables `u` for a given set of
`equations`, i.e. the waterheight plus the bathymetry.

`u` is a vector of the conserved variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.
"""
function waterheight_total end

"""
    waterheight(u, equations)

Return the waterheight of the conserved variables `u` for a given set of
`equations`, i.e. the waterheight above the bathymetry.

`u` is a vector of the conserved variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.
"""
function waterheight end

"""
    velocity(u, equations)

Return the velocity of the conserved variables `u` for a given set of
`equations`.

`u` is a vector of the conserved variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.
"""
function velocity end

"""
    momentum(u, equations)

Return the momentum of the conserved variables `u` for a given set of
`equations`, i.e. the waterheight times the velocity.

`u` is a vector of the conserved variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.
"""
@inline function momentum(u, equations::AbstractEquations)
    return waterheight(u, equations) * velocity(u, equations)
end

"""
    entropy(u, equations)

Return the chosen entropy of the conserved variables `u` for a given set of
`equations`.

`u` is a vector of the conserved variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.
"""
function entropy end

"""
    energy_total(u, equations)

Return the total energy of the conserved variables `u` for a given set of
`equations`.

`u` is a vector of the conserved variables at a single node, i.e., a vector
of the correct length `nvariables(equations)`.
"""
function energy_total end

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
            println("    variable " * string(variable), ": ", varnames(equations)[variable])
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

# Shallow water equations
abstract type AbstractShallowWaterEquations{NDIMS, NVARS} <: AbstractEquations{NDIMS, NVARS} end
include("shallow_water_1d.jl")

# BBM-BBM equations
abstract type AbstractBBMBBMEquations{NDIMS, NVARS} <: AbstractEquations{NDIMS, NVARS} end
include("bbm_bbm_1d.jl")
include("bbm_bbm_variable_bathymetry_1d.jl")
