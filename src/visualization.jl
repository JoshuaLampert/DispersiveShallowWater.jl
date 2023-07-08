struct PlotData{Ylim}
    semisol::Pair{<:Semidiscretization, <:ODESolution}
    plot_initial::Bool
    step::Integer
    # call this yli because ylim, ylimits etc. are already occupied by Plots.jl, but we want a vector
    yli::Ylim
end 

@recipe function f(plotdata::PlotData)
    @unpack semisol, plot_initial, step, yli = plotdata
    semi, sol = semisol
    nvars = nvariables(semi)

    # TODO: hardcode this for now, might need to adjuts this in the future
    nsubplots = 2
    if isnothing(yli)
        yli = fill(:auto, nsubplots)
    end
    if step == -1
        step = length(sol.t)
    end
    @assert length(yli)==nsubplots "The vector yli must be as long as there are variables"

    equations = semi.equations
    initial_condition = semi.initial_condition
    t = sol.t[step]

    if plot_initial == true
        u_exact = wrap_array(compute_coefficients(initial_condition, t, semi), semi)
    end

    data = wrap_array(sol.u[step], semi)
    bathy = zeros(nnodes(semi))
    for j in eachnode(semi)
        bathy[j] = bathymetry(view(data, :, j), equations)
    end

    names = (varnames(equations))
    plot_title --> "$(get_name(semi.equations)) at t = $(round(t, digits=5))"
    size --> (1200, 800)
    layout := nsubplots

    for i in 1:nsubplots
        if plot_initial == true
            @series begin
                subplot := i
                label := "initial $(names[i])"
                grid(semi), u_exact[i, :]
            end
        end

        # Don't plot bathymetry in separate subplot
        names[i] == "D" && continue
        @series begin
            subplot := i
            label := names[i]
            xguide := "x"
            yguide := names[i]
            title := names[i]
            ylim := yli[i]
            grid(semi), data[i, :]
        end
    end

    # Plot the bathymetry
    @series begin
        subplot := 1
        label := "bathymetry"
        xguide := "x"
        yguide := names[1]
        title := names[1]
        ylim := yli[1]
        color := :black
        grid(semi), bathy
    end
end

@recipe function f(semisol::Pair{<:Semidiscretization, <:ODESolution}; plot_initial = false,
                   step = -1, yli = nothing)
    PlotData(semisol, plot_initial, step, yli)
end

@recipe function f(semi::Semidiscretization, sol::ODESolution; plot_initial = false,
                   step = -1, yli = nothing)
    PlotData(semi => sol, plot_initial, step, yli)
end

# TODO: Only plot change in invariants for now, also plot errors?
@recipe function f(cb::DiscreteCallback{Condition, Affect!}) where {Condition, Affect! <: AnalysisCallback}
    t = tstops(cb)
    ints = integrals(cb)
    plot_title --> "change in invariants"
    for (name, integral) in pairs(ints)
        @series begin
            label := string(name)
            xguide := "t"
            yguide := "change in invariants"
            t, integral .- integral[1]
        end
    end
end
