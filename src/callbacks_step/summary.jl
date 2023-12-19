"""
    SummaryCallback(io::IO = stdout)

Create and return a callback that resets the timer at the beginning of
a simulation and prints the timer values at the end of the simulation.
"""
struct SummaryCallback
    io::IO

    function SummaryCallback(io::IO = stdout)
        function initialize(cb, u, t, integrator)
            initialize_summary_callback(cb, u, t, integrator)
        end
        # At the end of the simulation, the timer is printed
        function finalize(cb, u, t, integrator)
            finalize_summary_callback(cb, u, t, integrator)
        end
        summary_callback = new(io)
        # SummaryCallback is never called during the simulation
        condition = (u, t, integrator) -> false
        DiscreteCallback(condition, summary_callback,
                         save_positions = (false, false),
                         initialize = initialize,
                         finalize = finalize)
    end
end

function Base.show(io::IO, cb::DiscreteCallback{<:Any, <:SummaryCallback})
    @nospecialize cb # reduce precompilation time

    print(io, "SummaryCallback")
end

function initialize_summary_callback(cb::DiscreteCallback, u, t, integrator)
    reset_timer!(timer())
    return nothing
end

# the summary callback does nothing when called accidentally
(cb::SummaryCallback)(integrator) = u_modified!(integrator, false)

function finalize_summary_callback(cb::DiscreteCallback, u, t, integrator)
    io = cb.affect!.io
    TimerOutputs.complement!(timer())
    print_timer(io, timer(), title = "DispersiveSWE",
                allocations = true, linechars = :unicode, compact = false)
    println(io)
    return nothing
end
