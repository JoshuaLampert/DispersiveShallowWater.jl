"""
    examples_dir()

Return the directory where the example files provided with DispersiveShallowWater.jl are located. If DispersiveShallowWater is
installed as a regular package (with `]add DispersiveShallowWater`), these files are read-only and should *not* be
modified. To find out which files are available, use, e.g., `readdir`.

Copied from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).

# Examples
```@example
readdir(examples_dir())
```
"""
examples_dir() = pkgdir(DispersiveShallowWater, "examples")

"""
    get_examples()

Return a list of all examples that are provided by DispersiveShallowWater.jl. See also
[`examples_dir`](@ref) and [`default_example`](@ref).

Copied from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).
"""
function get_examples()
    examples = String[]
    for (root, dirs, files) in walkdir(examples_dir())
        for f in files
            if endswith(f, ".jl")
                push!(examples, joinpath(root, f))
            end
        end
    end

    return examples
end

"""
    default_example()

Return the path to an example that can be used to quickly see DispersiveShallowWater.jl in action.
See also [`examples_dir`](@ref) and [`get_examples`](@ref).

Copied from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).
"""
function default_example()
    joinpath(examples_dir(), "bbm_bbm_variable_bathymetry_1d",
             "bbm_bbm_variable_bathymetry_1d_basic.jl")
end

"""
    path_create_figures()

Return the path to the file that creates all figures used in the master thesis "Structure-preserving
Numerical Methods for Dispersive Shallow Water Model" (2023). Executing this julia script may take a
while.

# Examples
```@example
include(DispersiveShallowWater.path_create_figures())
```
"""
function path_create_figures()
    pkgdir(DispersiveShallowWater, "create_figures.jl")
end

# Note: We can't call the method below `DispersiveShallowWater.include` since that is created automatically
# inside `module DispersiveShallowWater` to `include` source files and evaluate them within the global scope
# of `DispersiveShallowWater`. However, users will want to evaluate in the global scope of `Main` or something
# similar to manage dependencies on their own.
"""
    trixi_include([mod::Module=Main,] example::AbstractString; kwargs...)

`include` the file `example` and evaluate its content in the global scope of module `mod`.
You can override specific assignments in `example` by supplying keyword arguments.
It's basic purpose is to make it easier to modify some parameters while running DispersiveShallowWater from the
REPL. Additionally, this is used in tests to reduce the computational burden for CI while still
providing examples with sensible default values for users.

Before replacing assignments in `example`, the keyword argument `maxiters` is inserted
into calls to `solve` and `DispersiveShallowWater.solve` with it's default value used in the SciML ecosystem
for ODEs, see the "Miscellaneous" section of the
[documentation](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/).

Copied from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).

# Examples

```jldoctest
julia> redirect_stdout(devnull) do
           trixi_include(@__MODULE__, joinpath(examples_dir(), "bbm_bbm_1d", "bbm_bbm_1d_basic.jl"),
                         tspan=(0.0, 0.1))
           sol.t[end]
       end
0.1
```
"""
function trixi_include(mod::Module, example::AbstractString; kwargs...)
    Base.include(ex -> replace_assignments(insert_maxiters(ex); kwargs...), mod, example)
end

trixi_include(example::AbstractString; kwargs...) = trixi_include(Main, example; kwargs...)

"""
    convergence_test([mod::Module=Main,] example::AbstractString, iterations; kwargs...)

Run `iterations` simulations using the setup given in `example` and compute
the experimental order of convergence (EOC) in the ``L^2`` and ``L^\\infty`` norm.
In each iteration, the resolution of the respective mesh will be doubled.
Additional keyword arguments `kwargs...` and the optional module `mod` are passed directly
to [`trixi_include`](@ref).

Adjusted from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).
"""
function convergence_test(mod::Module, example::AbstractString, iterations; kwargs...)
    @assert(iterations>1,
            "Number of iterations must be bigger than 1 for a convergence analysis")

    # Types of errors to be calculated
    errors = Dict(:l2 => Float64[], :linf => Float64[])

    initial_N = extract_initial_N(example, kwargs)

    # run simulations and extract errors
    for iter in 1:iterations
        println("Running convtest iteration ", iter, "/", iterations)

        trixi_include(mod, example; kwargs..., N = initial_N * 2^(iter - 1))

        l2_error, linf_error = mod.analysis_callback(mod.sol)

        # collect errors as one vector to reshape later
        append!(errors[:l2], l2_error)
        append!(errors[:linf], linf_error)

        println("\n\n")
        println("#"^100)
    end

    # Use raw error values to compute EOC
    analyze_convergence(errors, iterations, mod.semi, initial_N)
end

# Analyze convergence for any semidiscretization
# Note: this intermediate method is to allow dispatching on the semidiscretization
function analyze_convergence(errors, iterations, semi::Semidiscretization, initial_N)
    _, equations, _, _ = mesh_equations_solver_cache(semi)
    variablenames = varnames(prim2prim, equations)
    analyze_convergence(errors, iterations, variablenames, initial_N)
end

# This method is called with the collected error values to actually compute and print the EOC
function analyze_convergence(errors, iterations,
                             variablenames::Union{Tuple, AbstractArray}, initial_N)
    nvariables = length(variablenames)

    # Reshape errors to get a matrix where the i-th row represents the i-th iteration
    # and the j-th column represents the j-th variable
    errorsmatrix = Dict(kind => transpose(reshape(error, (nvariables, iterations)))
                        for (kind, error) in errors)

    # Calculate EOCs where the columns represent the variables
    # As dx halves in every iteration the denominator needs to be log(1/2)
    eocs = Dict(kind => log.(error[2:end, :] ./ error[1:(end - 1), :]) ./ log(1 / 2)
                for (kind, error) in errorsmatrix)

    eoc_mean_values = Dict{Symbol, Any}()
    eoc_mean_values[:variables] = variablenames

    for (kind, error) in errorsmatrix
        println(kind)

        for v in variablenames
            @printf("%-25s", v)
        end
        println("")

        for k in 1:nvariables
            @printf("%-5s", "N")
            @printf("%-10s", "error")
            @printf("%-10s", "EOC")
        end
        println("")

        # Print errors for the first iteration
        for k in 1:nvariables
            @printf("%-5d", initial_N)
            @printf("%-10.2e", error[1, k])
            @printf("%-10s", "-")
        end
        println("")

        # For the following iterations print errors and EOCs
        for j in 2:iterations
            for k in 1:nvariables
                @printf("%-5d", initial_N*2^(j - 1))
                @printf("%-10.2e", error[j, k])
                @printf("%-10.2f", eocs[kind][j - 1, k])
            end
            println("")
        end
        println("")

        # Print mean EOCs
        mean_values = zeros(nvariables)
        for v in 1:nvariables
            mean_values[v] = sum(eocs[kind][:, v]) ./ length(eocs[kind][:, v])
            @printf("%-10s", "mean")
            @printf("%-10.2f", mean_values[v])
        end
        eoc_mean_values[kind] = mean_values
        println("")
        println("-"^100)
    end

    return eoc_mean_values, errorsmatrix
end

function convergence_test(example::AbstractString, iterations; kwargs...)
    convergence_test(Main, example::AbstractString, iterations; kwargs...)
end

# Helper methods used in the functions defined above, also copied from Trixi.jl

# Apply the function `f` to `expr` and all sub-expressions recursively.
walkexpr(f, expr::Expr) = f(Expr(expr.head, (walkexpr(f, arg) for arg in expr.args)...))
walkexpr(f, x) = f(x)

# Insert the keyword argument `maxiters` into calls to `solve` and `DispersiveShallowWater.solve`
# with default value `10^5` if it is not already present.
function insert_maxiters(expr)
    maxiters_default = 10^5

    expr = walkexpr(expr) do x
        if x isa Expr
            is_plain_solve = x.head === Symbol("call") && x.args[1] === Symbol("solve")
            is_trixi_solve = (x.head === Symbol("call") && x.args[1] isa Expr &&
                              x.args[1].head === Symbol(".") &&
                              x.args[1].args[1] === Symbol("DispersiveShallowWater") &&
                              x.args[1].args[2] isa QuoteNode &&
                              x.args[1].args[2].value === Symbol("solve"))

            if is_plain_solve || is_trixi_solve
                # Do nothing if `maxiters` is already set as keyword argument...
                for arg in x.args
                    if arg isa Expr && arg.head === Symbol("kw") &&
                       arg.args[1] === Symbol("maxiters")
                        return x
                    end
                end

                # ...and insert it otherwise.
                push!(x.args, Expr(Symbol("kw"), Symbol("maxiters"), maxiters_default))
            end
        end
        return x
    end

    return expr
end

# Replace assignments to `key` in `expr` by `key = val` for all `(key,val)` in `kwargs`.
function replace_assignments(expr; kwargs...)
    # replace explicit and keyword assignments
    expr = walkexpr(expr) do x
        if x isa Expr
            for (key, val) in kwargs
                if (x.head === Symbol("=") || x.head === :kw) && x.args[1] === Symbol(key)
                    x.args[2] = :($val)
                    # dump(x)
                end
            end
        end
        return x
    end

    return expr
end

# find a (keyword or common) assignment to `destination` in `expr`
# and return the assigned value
function find_assignment(expr, destination)
    # declare result to be able to assign to it in the closure
    local result

    # find explicit and keyword assignments
    walkexpr(expr) do x
        if x isa Expr
            if (x.head === Symbol("=") || x.head === :kw) &&
               x.args[1] === Symbol(destination)
                result = x.args[2]
                # dump(x)
            end
        end
        return x
    end

    result
end

function extract_initial_N(example, kwargs)
    code = read(example, String)
    expr = Meta.parse("begin \n$code \nend")

    if haskey(kwargs, :N)
        return kwargs[:N]
    else
        # get N from the example
        N = find_assignment(expr, :N)
        return N
    end
end

"""
    @autoinfiltrate
    @autoinfiltrate condition::Bool

Invoke the `@infiltrate` macro of the package Infiltrator.jl to create a breakpoint for ad-hoc
interactive debugging in the REPL. If the optional argument `condition` is given, the breakpoint is
only enabled if `condition` evaluates to `true`.

As opposed to using `Infiltrator.@infiltrate` directly, this macro does not require Infiltrator.jl
to be added as a dependency to DispersiveShallowWater.jl. As a bonus, the macro will also attempt to load
the Infiltrator module if it has not yet been loaded manually.

Note: For this macro to work, the Infiltrator.jl package needs to be installed in your current Julia
environment stack.

See also: [Infiltrator.jl](https://github.com/JuliaDebug/Infiltrator.jl)

!!! warning "Internal use only"
    Please note that this macro is intended for internal use only. It is *not* part of the public
    API of DispersiveShallowWater.jl, and it thus can altered (or be removed) at any time without it being
    considered a breaking change.
"""
macro autoinfiltrate(condition = true)
    pkgid = Base.PkgId(Base.UUID("5903a43b-9cc3-4c30-8d17-598619ec4e9b"), "Infiltrator")
    if !haskey(Base.loaded_modules, pkgid)
        try
            Base.eval(Main, :(using Infiltrator))
        catch err
            @error "Cannot load Infiltrator.jl. Make sure it is included in your environment stack."
        end
    end
    i = get(Base.loaded_modules, pkgid, nothing)
    lnn = LineNumberNode(__source__.line, __source__.file)

    if i === nothing
        return Expr(:macrocall,
                    Symbol("@warn"),
                    lnn,
                    "Could not load Infiltrator.")
    end

    return Expr(:macrocall,
                Expr(:., i, QuoteNode(Symbol("@infiltrate"))),
                lnn,
                esc(condition))
end
