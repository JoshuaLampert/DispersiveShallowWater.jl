struct PlotData{Conversion}
    semisol::Pair{<:Semidiscretization, <:ODESolution}
    plot_initial::Bool
    plot_bathymetry::Bool
    conversion::Conversion
    step::Integer
end

struct PlotDataOverTime{RealT, Conversion}
    semisol::Pair{<:Semidiscretization, <:ODESolution}
    x::RealT
    conversion::Conversion
end

@recipe function f(plotdata::PlotData)
    @unpack semisol, plot_initial, plot_bathymetry, conversion, step = plotdata
    semi, sol = semisol
    equations = semi.equations
    names = varnames(conversion, equations)

    nvars = length(conversion(zeros(nvariables(semi)), equations))
    nsubplots = nvars
    # Bathymetry is not plotted in separate subplot
    if length(intersect(names, ["D", "b"])) > 0
        nsubplots -= 1
    end

    if step == -1
        step = length(sol.t)
    end

    initial_condition = semi.initial_condition
    t = sol.t[step]

    if plot_initial == true
        q_exact = compute_coefficients(initial_condition, t, semi)
        data_exact = zeros(nvars, nnodes(semi))
    end

    q = sol.u[step]
    data = zeros(nvars, nnodes(semi))
    if plot_bathymetry == true
        bathy = zeros(nnodes(semi))
    end
    for j in eachnode(semi)
        if plot_bathymetry == true
            bathy[j] = bathymetry(get_node_vars(q, equations, j), equations)
        end
        if plot_initial == true
            data_exact[:, j] .= conversion(get_node_vars(q_exact, equations, j), equations)
        end
        data[:, j] .= conversion(get_node_vars(q, equations, j), equations)
    end

    plot_title --> "$(get_name(equations)) at t = $(round(t, digits=5))"
    layout --> nsubplots

    for i in 1:nsubplots
        # Don't plot bathymetry in separate subplot
        names[i] in ["D", "b"] && continue

        if plot_initial == true
            @series begin
                subplot --> i
                linestyle := :solid
                label := "initial $(names[i])"
                grid(semi), data_exact[i, :]
            end
        end

        @series begin
            subplot --> i
            label --> names[i]
            xguide --> "x"
            yguide --> names[i]
            title --> names[i]
            grid(semi), data[i, :]
        end
    end

    # Plot the bathymetry
    if plot_bathymetry == true
        @series begin
            subplot --> 1
            linestyle := :solid
            label := "bathymetry"
            xguide --> "x"
            yguide --> names[1]
            title --> names[1]
            color := :black
            grid(semi), bathy
        end
    end
end

@recipe function f(plotdata_over_time::PlotDataOverTime)
    @unpack semisol, x, conversion = plotdata_over_time
    semi, sol = semisol
    equations = semi.equations
    names = varnames(conversion, equations)

    nvars = length(conversion(zeros(nvariables(semi)), equations))
    nsubplots = nvars
    # Bathymetry is not plotted in separate subplot
    if length(intersect(names, ["D", "b"])) > 0
        nsubplots -= 1
    end

    solution = zeros(nvariables(semi), length(sol.t))
    data = zeros(nvars, length(sol.t))
    for v in eachvariable(semi)
        for k in 1:length(sol.t)
            # Allow that the spatial value `x` is not on the grid. Thus, interpolate the given values to the provided `x`
            # with a linear spline.
            solution[v, k] = linear_interpolation(grid(semi),
                                                  sol.u[k].x[v])(x)
        end
    end

    for k in 1:length(sol.t)
        data[:, k] .= conversion(view(solution, :, k), equations)
    end

    names = varnames(conversion, equations)
    plot_title -->
    "$(get_name(semi.equations)) at x = $(round(x, digits=5))"
    layout --> nsubplots

    for i in 1:nsubplots
        # Don't plot bathymetry in separate subplot
        names[i] in ["D", "b"] && continue
        @series begin
            subplot --> i
            label --> names[i]
            xguide --> "t"
            yguide --> names[i]
            title --> names[i]
            sol.t, data[i, :]
        end
    end
end

@recipe function f(semisol::Pair{<:Semidiscretization, <:ODESolution}; plot_initial = false,
                   plot_bathymetry = true, conversion = prim2phys, step = -1)
    PlotData(semisol, plot_initial, plot_bathymetry, conversion, step)
end

@recipe function f(semi::Semidiscretization, sol::ODESolution; plot_initial = false,
                   plot_bathymetry = true, conversion = prim2phys, step = -1)
    PlotData(semi => sol, plot_initial, plot_bathymetry, conversion, step)
end

@recipe function f(semisol::Pair{<:Semidiscretization, <:ODESolution}, x_value;
                   conversion = prim2phys)
    PlotDataOverTime(semisol, x_value, conversion)
end

@recipe function f(semi::Semidiscretization, sol::ODESolution, x_value;
                   conversion = prim2phys)
    PlotDataOverTime(semi => sol, x_value, conversion)
end

function pretty(name)
    if name == :waterheight_total
        return "∫η"
    elseif name == :velocity
        return "∫v"
    elseif name in [:discharge, :momentum]
        return "∫P"
    elseif name == :entropy
        return "∫U"
    elseif name == :entropy_modified
        return "∫U_mod"
    elseif name == :l2_error
        return "L² error"
    elseif name == :linf_error
        return "L∞ error"
    elseif name == :conservation_error
        return "∫|q_q₀|"
    elseif name == :lake_at_rest_error
        return "∫|η-η₀|"
    else
        return string(name)
    end
end

@recipe function f(cb::DiscreteCallback{Condition, Affect!}; what = (:integrals,),
                   label_extension = "", start_from = 1,
                   exclude = []) where {Condition, Affect! <: AnalysisCallback}
    t = tstops(cb)
    @assert length(t)>start_from "The keyword argument `start_from` needs to be smaller than the number of timesteps: $(length(t))"
    subplot = 1
    layout --> length(what)
    if :integrals in what
        ints = integrals(cb)

        for (i, (name, integral)) in enumerate(pairs(ints))
            name in exclude && continue
            @series begin
                subplot --> subplot
                label := pretty(name) * " " * label_extension
                title --> "change of invariants"
                xguide --> "t"
                yguide --> "change of invariants"
                t[start_from:end], (integral .- integral[1])[start_from:end]
            end
        end
        subplot += 1
    end
    if :errors in what
        errs = errors(cb)
        for (i, (name, err)) in enumerate(pairs(errs))
            name in exclude && continue
            @series begin
                subplot --> subplot
                label --> pretty(name) * " " * label_extension
                title --> "errors"
                xguide --> "t"
                yguide --> "sum of errors"
                t[start_from:end], dropdims(sum(err, dims = 1), dims = 1)[start_from:end]
            end
        end
    end
end
