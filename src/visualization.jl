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
    equations = semi.equations
    names = varnames(conversion, equations)

    nvars = length(conversion(zeros(nvariables(semi)), equations))
    nsubplots = nvars
    # Bathymetry is not plotted in separate subplot
    if length(intersect(names, ["D", "b"])) > 0
        nsubplots -= 1
    end

    if isnothing(yli)
        yli = fill(:auto, nsubplots)
    end
    if step == -1
        step = length(sol.t)
    end
    @assert length(yli)==nsubplots "The vector yli must be as long as there are subplots"

    initial_condition = semi.initial_condition
    t = sol.t[step]

    if plot_initial == true
        q_exact = wrap_array(compute_coefficients(initial_condition, t, semi), semi)
    end

    q = wrap_array(sol.u[step], semi)
    data = zeros(nvars, nnodes(semi))
    if plot_bathymetry == true
        bathy = zeros(nnodes(semi))
    end
    for j in eachnode(semi)
        if plot_bathymetry == true
            bathy[j] = bathymetry(view(q, :, j), equations)
        end
        if plot_initial == true
            q_exact[:, j] .= conversion(view(q_exact, :, j), equations)
        end
        data[:, j] .= conversion(view(q, :, j), equations)
    end

    plot_title --> "$(get_name(semi.equations)) at t = $(round(t, digits=5))"
    size --> (1200, 800)
    layout := nsubplots

    for i in 1:nsubplots
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

    # Plot the bathymetry
    if plot_bathymetry == true
        @series begin
            subplot := 1
            linestyle := :solid
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
    equations = semi.equations
    names = varnames(conversion, equations)

    nvars = length(conversion(zeros(nvariables(semi)), equations))
    nsubplots = nvars
    # Bathymetry is not plotted in separate subplot
    if length(intersect(names, ["D", "b"])) > 0
        nsubplots -= 1
    end

    if isnothing(yli)
        yli = fill(:auto, nsubplots)
    end

    index = argmin(abs.(grid(semi) .- x))
    solution = zeros(nvariables(semi), length(sol.t))
    data = zeros(nvars, length(sol.t))
    for i in 1:nvariables(semi)
        for k in 1:length(sol.t)
            solution[i, k] = wrap_array(sol.u[k], semi)[i, index]
        end
    end

    for k in 1:length(sol.t)
        data[:, k] .= conversion(view(solution, :, k), equations)
    end

    names = varnames(conversion, equations)
    plot_title -->
    "$(get_name(semi.equations)) at x = $(round(grid(semi)[index], digits=5))"
    size --> (1200, 800)
    layout := nsubplots

    for i in 1:nsubplots
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

@recipe function f(cb::DiscreteCallback{Condition, Affect!}; what = (:integrals,),
                   label_extension = "") where {Condition, Affect! <: AnalysisCallback}
    t = tstops(cb)
    subplot = 1
    layout := length(what)
    if :integrals in what
        ints = integrals(cb)
        plot_title --> "change of invariants"
        for (name, integral) in pairs(ints)
            @series begin
                subplot := subplot
                label := string(name) * " " * label_extension
                xguide --> "t"
                yguide --> "change of invariants"
                t, integral .- integral[1]
            end
        end
        subplot += 1
    end
    if :errors in what
        errs = errors(cb)
        plot_title --> "errors"
        for (name, err) in pairs(errs)
            @series begin
                subplot := subplot
                label := string(name) * " " * label_extension
                xguide --> "t"
                yguide --> "sum of errors"
                t, dropdims(sum(err, dims = 1), dims = 1)
            end
        end
    end
end
