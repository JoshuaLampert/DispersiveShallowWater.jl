struct PlotData{Ylim, Conversion}
    semisol::Pair{<:Semidiscretization, <:ODESolution}
    plot_initial::Bool
    plot_bathymetry::Bool
    conversion::Conversion
    step::Integer
    # call this yli because ylim, ylimits etc. are already occupied by Plots.jl, but we want a vector
    yli::Ylim
end

struct PlotDataOverTime{RealT, Ylim, Conversion}
    semisol::Pair{<:Semidiscretization, <:ODESolution}
    x::RealT
    conversion::Conversion
    # call this yli because ylim, ylimits etc. are already occupied by Plots.jl, but we want a vector
    yli::Ylim
end

@recipe function f(plotdata::PlotData)
    @unpack semisol, plot_initial, plot_bathymetry, conversion, step, yli = plotdata
    semi, sol = semisol
    nvars = nvariables(semi)

    # TODO: hardcode this for now, might need to adjust this in the future
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
        q_exact = wrap_array(compute_coefficients(initial_condition, t, semi), semi)
    end

    q = wrap_array(sol.u[step], semi)
    data = similar(q)
    if conversion == prim2prim && plot_bathymetry == true
        bathy = zeros(nnodes(semi))
    end
    for j in eachnode(semi)
        if conversion == prim2prim && plot_bathymetry == true
            bathy[j] = bathymetry(view(q, :, j), equations)
        end
        if plot_initial == true
            q_exact[:, j] = conversion(view(q_exact, :, j), equations)
        end
        data[:, j] = conversion(view(q, :, j), equations)
    end

    names = varnames(conversion, equations)
    plot_title --> "$(get_name(semi.equations)) at t = $(round(t, digits=5))"
    size --> (1200, 800)
    layout := nsubplots

    for i in 1:nvars
        # Don't plot bathymetry in separate subplot
        names[i] in ["D", "b"] && continue

        if plot_initial == true
            @series begin
                subplot := i
                label := "initial $(names[i])"
                grid(semi), q_exact[i, :]
            end
        end

        @series begin
            subplot := i
            label --> names[i]
            xguide := "x"
            yguide := names[i]
            title := names[i]
            ylim := yli[i]
            grid(semi), data[i, :]
        end
    end

    # Plot the bathymetry if primitive variables are plotted
    if conversion == prim2prim && plot_bathymetry == true
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
end

@recipe function f(plotdata_over_time::PlotDataOverTime)
    @unpack semisol, x, conversion, yli = plotdata_over_time
    semi, sol = semisol
    nvars = nvariables(semi)

    # TODO: hardcode this for now, might need to adjust this in the future
    nsubplots = 2
    if isnothing(yli)
        yli = fill(:auto, nsubplots)
    end
    equations = semi.equations

    index = argmin(abs.(grid(semi) .- x))
    data = zeros(nvars, length(sol.t))
    for i in 1:nvars
        for k in 1:length(sol.t)
            data[i, k] = wrap_array(sol.u[k], semi)[i, index]
        end
    end

    for k in 1:length(sol.t)
        data[:, k] = conversion(view(data, :, k), equations)
    end

    names = varnames(conversion, equations)
    plot_title -->
    "$(get_name(semi.equations)) at x = $(round(grid(semi)[index], digits=5))"
    size --> (1200, 800)
    layout := nsubplots

    for i in 1:nvars
        # Don't plot bathymetry in separate subplot
        names[i] in ["D", "b"] && continue
        @series begin
            subplot := i
            label := names[i]
            xguide := "t"
            yguide := names[i]
            title := names[i]
            ylim := yli[i]
            sol.t, data[i, :]
        end
    end
end

@recipe function f(semisol::Pair{<:Semidiscretization, <:ODESolution}; plot_initial = false,
                   plot_bathymetry = true, conversion = prim2prim, step = -1, yli = nothing)
    PlotData(semisol, plot_initial, plot_bathymetry, conversion, step, yli)
end

@recipe function f(semi::Semidiscretization, sol::ODESolution; plot_initial = false,
                   plot_bathymetry = true, conversion = prim2prim, step = -1, yli = nothing)
    PlotData(semi => sol, plot_initial, plot_bathymetry, conversion, step, yli)
end

@recipe function f(semisol::Pair{<:Semidiscretization, <:ODESolution}, x_value;
                   conversion = prim2prim, yli = nothing)
    PlotDataOverTime(semisol, x_value, conversion, yli)
end

@recipe function f(semi::Semidiscretization, sol::ODESolution, x_value;
                   conversion = prim2prim, yli = nothing)
    PlotDataOverTime(semi => sol, x_value, conversion, yli)
end

# TODO: Only plot change in invariants for now, also plot errors?
@recipe function f(cb::DiscreteCallback{Condition, Affect!}) where {Condition,
                                                                    Affect! <:
                                                                    AnalysisCallback}
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
