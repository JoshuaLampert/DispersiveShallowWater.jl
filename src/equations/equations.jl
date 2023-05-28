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

@inline Base.ndims(::AbstractEquations{NDIMS}) where NDIMS = NDIMS


# BBM-BBM equations
abstract type AbstractBBMBBMEquations{NDIMS, NVARS} <: AbstractEquations{NDIMS, NVARS} end
include("bbm_bbm_1d.jl")
