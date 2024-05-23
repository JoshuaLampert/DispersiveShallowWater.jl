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

function convergence_test(example::AbstractString, iterations_or_Ns; kwargs...)
    convergence_test(Main, example::AbstractString, iterations_or_Ns; kwargs...)
end

"""
    convergence_test([mod::Module=Main,] example::AbstractString, iterations; kwargs...)
    convergence_test([mod::Module=Main,] example::AbstractString, Ns::AbstractVector; kwargs...)

Run multiple simulations using the setup given in `example` and compute
the experimental order of convergence (EOC) in the ``L^2`` and ``L^\\infty`` norm.
If `iterations` is passed as integer, in each iteration, the resolution of the respective mesh
will be doubled. If `Ns` is passed as vector, the simulations will be run for each value of `Ns`.
Additional keyword arguments `kwargs...` and the optional module `mod` are passed directly
to [`trixi_include`](@ref).

Adjusted from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).
"""
function convergence_test(mod::Module, example::AbstractString, iterations; kwargs...)
    @assert(iterations>1,
            "Number of iterations must be bigger than 1 for a convergence analysis")

    initial_N = extract_initial_N(example, kwargs)
    Ns = initial_N * 2 .^ (0:(iterations - 1))
    convergence_test(mod, example, Ns; kwargs...)
end

function convergence_test(mod::Module, example::AbstractString, Ns::AbstractVector;
                          kwargs...)
    # Types of errors to be calculated
    errors = Dict(:l2 => Float64[], :linf => Float64[])

    sort!(Ns)
    iterations = length(Ns)
    # run simulations and extract errors
    for iter in 1:iterations
        println("Running convtest iteration ", iter, "/", iterations)

        trixi_include(mod, example; kwargs..., N = Ns[iter])

        l2_error, linf_error = mod.analysis_callback(mod.sol)

        # collect errors as one vector to reshape later
        append!(errors[:l2], l2_error)
        append!(errors[:linf], linf_error)

        println("\n\n")
        println("#"^100)
    end

    # Use raw error values to compute EOC
    analyze_convergence(errors, iterations, mod.semi, Ns)
end

# Analyze convergence for any semidiscretization
# Note: this intermediate method is to allow dispatching on the semidiscretization
function analyze_convergence(errors, iterations, semi::Semidiscretization, Ns)
    _, equations, _, _ = mesh_equations_solver_cache(semi)
    variablenames = varnames(prim2prim, equations)
    analyze_convergence(errors, iterations, variablenames, Ns)
end

# This method is called with the collected error values to actually compute and print the EOC
function analyze_convergence(errors, iterations,
                             variablenames::Union{Tuple, AbstractArray}, Ns)
    nvariables = length(variablenames)

    # Reshape errors to get a matrix where the i-th row represents the i-th iteration
    # and the j-th column represents the j-th variable
    errorsmatrix = Dict(kind => transpose(reshape(error, (nvariables, iterations)))
                        for (kind, error) in errors)

    # Calculate EOCs where the columns represent the variables
    eocs = Dict(kind => log.(error[2:end, :] ./ error[1:(end - 1), :]) ./
                        log(Ns[1:(end - 1)] ./ Ns[2:end])
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
            @printf("%-5d", Ns[1])
            @printf("%-10.2e", error[1, k])
            @printf("%-10s", "-")
        end
        println("")

        # For the following iterations print errors and EOCs
        for j in 2:iterations
            for k in 1:nvariables
                @printf("%-5d", Ns[j])
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
            @printf("%-15s", "mean")
            @printf("%-10.2f", mean_values[v])
        end
        eoc_mean_values[kind] = mean_values
        println("")
        println("-"^100)
    end

    return eoc_mean_values, errorsmatrix
end

function extract_initial_N(example, kwargs)
    code = read(example, String)
    expr = Meta.parse("begin \n$code \nend")

    if haskey(kwargs, :N)
        return kwargs[:N]
    else
        # get N from the example
        N = TrixiBase.find_assignment(expr, :N)
        return N
    end
end

# Store main timer for global timing of functions
const main_timer = TimerOutput()

# Always call timer() to hide implementation details
timer() = main_timer

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
