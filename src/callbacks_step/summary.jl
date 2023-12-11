"""
    SummaryCallback(io::IO = stdout)

Create and return a callback that prints a human-readable summary of the simulation setup at the
beginning of a simulation and then resets the timer. When the returned callback is executed
directly, the current timer values are shown.
"""
struct SummaryCallback
    io::IO

    function SummaryCallback(io::IO = stdout)
        function initialize(cb, u, t, integrator)
            initialize_summary_callback(cb, u, t, integrator)
        end
        summary_callback = new(io)
        # SummaryCallback is called at end of simulation
        condition = (u, t, integrator) -> isfinished(integrator)
        DiscreteCallback(condition, summary_callback,
                         save_positions = (false, false),
                         initialize = initialize)
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

function (cb::SummaryCallback)(integrator)
    u_modified!(integrator, false)
    cb()
end

function (summary_callback::SummaryCallback)()
    io = summary_callback.io
    TimerOutputs.complement!(timer())
    print_timer(io, timer(), title = "DispersiveSWE",
                allocations = true, linechars = :unicode, compact = false)
    println(io)
    return nothing
end

# Allow calling the callback explicitly
function (cb::DiscreteCallback{Condition, Affect!})() where {Condition,
                                                             Affect! <: SummaryCallback}
    cb.affect!()
end
