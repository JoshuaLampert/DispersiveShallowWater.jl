"""
    AnalysisCallback(semi; interval=0,
                           extra_analysis_errors=Symbol[],
                           extra_analysis_integrals=(),
                           io=stdout)

Analyze a numerical solution every `interval` time steps.
The L2- and the L∞-norm for each component are computed by default.
Additional errors can be computed, e.g. by passing `extra_analysis_errors = (:conservation_error,)`.

Further scalar functions `func` in `extra_analysis_integrals` are applied to the numerical
solution and integrated over the computational domain. Some examples for this are
[`entropy`](@ref), and [`energy_total`](@ref).
You can also write your own function with the same signature as the examples listed above and
pass it via `extra_analysis_integrals`.
The computed errors and intergrals are saved for each timestep and can be obtained by calling
[`errors`](@ref) and [`integrals`](@ref).

During the Simulation, the `AnalysisCallback` will print information to `io`.
"""
mutable struct AnalysisCallback{T, AnalysisIntegrals, InitialStateIntegrals}
    start_time::Float64
    start_gc_time::Float64
    interval::Int
    analysis_errors::Vector{Symbol}
    analysis_integrals::AnalysisIntegrals
    initial_state_integrals::InitialStateIntegrals
    tstops::Vector{Float64}
    errors::Vector{Matrix{T}}
    integrals::Vector{Vector{T}}
    io::IO
end

function Base.show(io::IO, cb::DiscreteCallback{<:Any, <:AnalysisCallback})
    @nospecialize cb # reduce precompilation time

    analysis_callback = cb.affect!
    @unpack interval = analysis_callback
    print(io, "AnalysisCallback(interval=", interval, ")")
end

function Base.show(io::IO, ::MIME"text/plain",
                   cb::DiscreteCallback{<:Any, <:AnalysisCallback})
    @nospecialize cb # reduce precompilation time

    if get(io, :compact, false)
        show(io, cb)
    else
        analysis_callback = cb.affect!
        @unpack interval, analysis_errors, analysis_integrals = analysis_callback

        println(io, "AnalysisCallback")
        println(io, "    interval: ", interval)
        for (idx, error) in enumerate(analysis_errors)
            println(io, "    error ", idx, ": ", error)
        end
        for (idx, integral) in enumerate(analysis_integrals)
            println(io, "    integral ", idx, ": ", integral)
        end
    end
end

function AnalysisCallback(semi::Semidiscretization; kwargs...)
    mesh, equations, solver = mesh_equations_solver(semi)
    AnalysisCallback(mesh, equations, solver; kwargs...)
end

function AnalysisCallback(mesh, equations::AbstractEquations, solver;
                          interval = 0,
                          extra_analysis_errors = Symbol[],
                          analysis_errors = union(default_analysis_errors(equations),
                                                  extra_analysis_errors),
                          extra_analysis_integrals = (),
                          analysis_integrals = union(default_analysis_integrals(equations),
                                                     extra_analysis_integrals),
                          io = stdout)
    # Decide when the callback is activated.
    # With error-based step size control, some steps can be rejected. Thus,
    #   `integrator.iter >= integrator.stats.naccept`
    #    (total #steps)       (#accepted steps)
    # We need to check the number of accepted steps since callbacks are not
    # activated after a rejected step.
    condition = (u, t, integrator) -> interval > 0 &&
        ((integrator.stats.naccept % interval == 0 &&
          !(integrator.stats.naccept == 0 && integrator.iter > 0)) ||
         isfinished(integrator))
    for extra_analysis_error in extra_analysis_errors
        if extra_analysis_error != :conservation_error
            @warn "Extra analysis error $extra_analysis_error is not supported and will be ignored."
        end
    end
    analysis_callback = AnalysisCallback(0.0, 0.0, interval, analysis_errors,
                                         Tuple(analysis_integrals),
                                         SVector(ntuple(_ -> zero(real(solver)),
                                                        Val(nvariables(equations)))),
                                         Vector{Float64}(),
                                         Vector{Matrix{real(solver)}}(),
                                         Vector{Vector{real(solver)}}(),
                                         io)

    DiscreteCallback(condition, analysis_callback,
                     save_positions = (false, false),
                     initialize = initialize!)
end

"""
    tstops(analysis_callback)

Return the time values that correspond to the saved values of the [`errors`](@ref) and [`integrals`](@ref).
"""
function tstops(cb::DiscreteCallback{
                                     Condition,
                                     Affect!
                                     }) where {Condition,
                                               Affect! <:
                                               AnalysisCallback}
    analysis_callback = cb.affect!
    return analysis_callback.tstops
end

"""
    errors(analysis_callback)
    
Return the computed errors for each timestep as a named tuple.
The shape of each entry is (nvariables, ntimesteps).
"""
function errors(cb::DiscreteCallback{
                                     Condition,
                                     Affect!
                                     }) where {Condition,
                                               Affect! <:
                                               AnalysisCallback}
    analysis_callback = cb.affect!
    names = collect(analysis_callback.analysis_errors)
    # "transpose" vector of matrices, first write it as 3d array and then convert it back to vector of matrices
    errors_array = reshape(reduce(hcat, analysis_callback.errors),
                           size(analysis_callback.errors[1])..., :)
    errors_vector = [errors_array[i, :, :] for i in 1:size(errors_array)[1]]
    return (; zip(names, errors_vector)...)
end

"""
    integrals(analysis_callback)

Return the computed integrals for each timestep as a named tuple.
"""
function integrals(cb::DiscreteCallback{
                                        Condition,
                                        Affect!
                                        }) where {Condition,
                                                  Affect! <:
                                                  AnalysisCallback}
    analysis_callback = cb.affect!
    names = collect(Symbol.(nameof.(analysis_callback.analysis_integrals)))
    # "transpose" vector of vector
    integrals = collect(eachrow(reduce(hcat, analysis_callback.integrals)))
    return (; zip(names, integrals)...)
end

function initialize!(cb::DiscreteCallback{Condition, Affect!}, u_ode, t,
                     integrator) where {Condition, Affect! <: AnalysisCallback}
    semi = integrator.p
    initial_state_integrals = integrate(u_ode, semi)

    analysis_callback = cb.affect!
    analysis_callback.initial_state_integrals = initial_state_integrals

    # Record current time using a high-resolution clock
    analysis_callback.start_time = time_ns()

    # Record total time spent in garbage collection so far using a high-resolution clock
    # Note: For details see the actual callback function below
    analysis_callback.start_gc_time = Base.gc_time_ns()

    analysis_callback(integrator)
    return nothing
end

function (analysis_callback::AnalysisCallback)(integrator)
    semi = integrator.p
    mesh, equations, solver = mesh_equations_solver(semi)

    l2_error, linf_error = analysis_callback(analysis_callback.io, integrator.u, integrator, semi)

    # avoid re-evaluating possible FSAL stages
    u_modified!(integrator, false)

    # Return errors for EOC analysis
    return l2_error, linf_error
end

# This method is just called internally from `(analysis_callback::AnalysisCallback)(integrator)`
# and serves as a function barrier. Additionally, it makes the code easier to profile and optimize.
function (analysis_callback::AnalysisCallback)(io, u_ode, integrator, semi)
    _, equations, solver = mesh_equations_solver(semi)
    @unpack analysis_errors, analysis_integrals, tstops, errors, integrals = analysis_callback
    @unpack t, dt = integrator
    push!(tstops, t)
    t_final = integrator.sol.prob.tspan[2]
    iter = integrator.stats.naccept

    # Compute the total runtime since the analysis callback has been initialized, in seconds
    runtime_absolute = 1.0e-9 * (time_ns() - analysis_callback.start_time)

    # Compute the total time spent in garbage collection since the analysis callback has been
    # initialized, in seconds
    # Note: `Base.gc_time_ns()` is not part of the public Julia API but has been available at least
    #        since Julia 1.6. Should this function be removed without replacement in a future Julia
    #        release, just delete this analysis quantity from the callback.
    # Source: https://github.com/JuliaLang/julia/blob/b540315cb4bd91e6f3a3e4ab8129a58556947628/base/timing.jl#L83-L84
    gc_time_absolute = 1.0e-9 * (Base.gc_time_ns() - analysis_callback.start_gc_time)

    # Compute the percentage of total time that was spent in garbage collection
    gc_time_percentage = gc_time_absolute / runtime_absolute * 100

    # Obtain the current memory usage of the Julia garbage collector, in MiB, i.e., the total size of
    # objects in memory that have been allocated by the JIT compiler or the user code.
    # Note: `Base.gc_live_bytes()` is not part of the public Julia API but has been available at least
    #        since Julia 1.6. Should this function be removed without replacement in a future Julia
    #        release, just delete this analysis quantity from the callback.
    # Source: https://github.com/JuliaLang/julia/blob/b540315cb4bd91e6f3a3e4ab8129a58556947628/base/timing.jl#L86-L97
    memory_use = Base.gc_live_bytes() / 2^20 # bytes -> MiB

    println(io)
    println(io, "─"^100)
    println(io, "Simulation running '", get_name(equations), "' with '", semi.initial_condition,
            "'")
    println(io, "─"^100)
    println(io, " #timesteps:     " * @sprintf("% 14d", iter) *
            "               " *
            " run time:       " * @sprintf("%10.8e s", runtime_absolute))
    println(io, " Δt:             " * @sprintf("%10.8e", dt) *
            "               " *
            " └── GC time:    " *
            @sprintf("%10.8e s (%5.3f%%)", gc_time_absolute, gc_time_percentage))
    println(io, " sim. time:      " * @sprintf("%10.8e (%5.3f%%)", t, t / t_final*100))
    println(io, " #DOF:           " * @sprintf("% 14d", nnodes(semi)) *
            "               " *
            " alloc'd memory: " * @sprintf("%14.3f MiB", memory_use))
    println(io)

    print(io, " Variable:    ")
    for v in eachvariable(equations)
        @printf(io, "   %-14s", varnames(prim2prim, equations)[v])
    end
    println(io)

    # Calculate L2/Linf errors, which are also returned
    l2_error, linf_error = calc_error_norms(u_ode, t, semi)
    current_errors = zeros(real(semi), (length(analysis_errors), nvariables(equations)))
    current_errors[1, :] = l2_error
    current_errors[2, :] = linf_error
    print(io, " L2 error:    ")
    for v in eachvariable(equations)
        @printf(io, "  % 10.8e", l2_error[v])
    end
    println(io)

    print(io, " Linf error:  ")
    for v in eachvariable(equations)
        @printf(io, "  % 10.8e", linf_error[v])
    end
    println(io)

    # Conservation error
    if :conservation_error in analysis_errors
        @unpack initial_state_integrals = analysis_callback
        state_integrals = integrate(u_ode, semi)
        current_errors[3, :] = abs.(state_integrals - initial_state_integrals)
        print(io, " |∫q - ∫q₀|:  ")
        for v in eachvariable(equations)
            @printf(io, "  % 10.8e", current_errors[3, v])
        end
        println(io)
    end
    push!(errors, current_errors)

    # additional integrals
    if length(analysis_integrals) > 0
        println(io)
        println(io, " Integrals:    ")
    end
    current_integrals = zeros(real(semi), length(analysis_integrals))
    analyze_integrals!(io, current_integrals, 1, analysis_integrals, u_ode, t, semi)
    push!(integrals, current_integrals)

    println(io, "─"^100)
    return l2_error, linf_error
end

# Iterate over tuples of analysis integrals in a type-stable way using "lispy tuple programming".
function analyze_integrals!(io, current_integrals, i, analysis_integrals::NTuple{N, Any}, u_ode,
                            t, semi) where {N}

    # Extract the first analysis integral and process it; keep the remaining to be processed later
    quantity = first(analysis_integrals)
    remaining_quantities = Base.tail(analysis_integrals)

    res = analyze(quantity, u_ode, t, semi)
    current_integrals[i] = res
    @printf(io, " %-12s:", pretty_form_utf(quantity))
    @printf(io, "  % 10.8e", res)
    println(io)

    # Recursively call this method with the unprocessed integrals
    analyze_integrals!(io, current_integrals, i + 1, remaining_quantities, u_ode, t, semi)
    return nothing
end

# terminate the type-stable iteration over tuples
function analyze_integrals!(io, current_integrals, i, analysis_integrals::Tuple{}, u_ode, t,
                            semi)
    nothing
end

# used for error checks and EOC analysis
function (cb::DiscreteCallback{Condition, Affect!})(sol) where {Condition,
                                                                Affect! <:
                                                                AnalysisCallback}
    analysis_callback = cb.affect!
    semi = sol.prob.p

    l2_error, linf_error = calc_error_norms(sol.u[end], sol.t[end], semi)

    return (; l2 = l2_error, linf = linf_error)
end

function analyze(quantity, u_ode, t, semi::Semidiscretization)
    integrate_quantity(u_ode -> quantity(u_ode, semi.equations), u_ode, semi)
end

# modified entropy from Svärd-Kalisch equations need to take the whole vector `u` for every point in space
function analyze(quantity::Union{typeof(energy_total_modified), typeof(entropy_modified)},
                 u_ode, t, semi::Semidiscretization)
    integrate_quantity(quantity, u_ode, semi)
end

pretty_form_utf(::typeof(waterheight_total)) = "∫η"
pretty_form_utf(::typeof(velocity)) = "∫v"
pretty_form_utf(::typeof(momentum)) = "∫P"
pretty_form_utf(::typeof(discharge)) = "∫P"
pretty_form_utf(::typeof(entropy)) = "∫U"
pretty_form_utf(::typeof(energy_total)) = "∫e_total"
pretty_form_utf(::typeof(entropy_modified)) = "∫U_modified"
pretty_form_utf(::typeof(energy_total_modified)) = "∫e_modified"
pretty_form_utf(::typeof(lake_at_rest_error)) = "∫|η-η₀|"
