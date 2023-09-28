using DispersiveShallowWater
using SummationByPartsOperators: legendre_derivative_operator,
                                 UniformPeriodicMesh1D,
                                 couple_discontinuously,
                                 couple_continuously
using SparseArrays: sparse
using Plots
using LaTeXStrings

const OUT = "out/"
ispath(OUT) || mkpath(OUT)
const EXAMPLES_DIR_BBMBBM = "bbm_bbm_1d"
const EXAMPLES_DIR_BBMBBM_VARIABLE = "bbm_bbm_variable_bathymetry_1d"
const EXAMPLES_DIR_SVAERD_KALISCH = "svaerd_kalisch_1d"

function run_example(path; kwargs...)
    trixi_include(joinpath(examples_dir(), path); kwargs...)
end

# Chapter 2
# Plot of bathymetry and waterheight
function fig_1()
    L = 1.0
    n = 100
    x = LinRange(0.0, L, n)
    fontsize = 20

    # just pick some function for b and eta that look nice
    H = 1.012

    b(x) = x * cos.(3 * pi * x) + H
    plot(x, b, color = :gray, fill = (0, 0.8, :gray), fillstyle = :/, linewidth = 3,
         legend = nothing, ticks = nothing, border = :none)

    eta(x) = x / (x^2 + 1) * sin(2 * pi * x) + H + 1.5
    plot!(x, eta, color = :blue, fill = (b.(x), 0.4, :blue), linewidth = 3)

    x1 = 0.2
    plot!([x1, x1], [b(x1), eta(x1)], line = (Plots.Arrow(:open, :both, 2.5, 2.0), :black),
          annotation = (x1 - 0.08, (eta(x1) + b(x1)) / 2, text(L"h(t, x)", fontsize)),
          linewidth = 2)
    x2 = 0.4
    plot!([x2, x2], [0.0, b(x2)], line = (Plots.Arrow(:open, :both, 2.5, 2.0), :black),
          annotation = (x2 + 0.06, b(x2) / 2, text(L"b(x)", fontsize)), linewidth = 2)
    x3 = 0.8
    plot!([x3, x3], [0.0, eta(x3)], line = (Plots.Arrow(:open, :both, 2.5, 2.0), :black),
          annotation = (x3 - 0.08, eta(x3) / 2, text(L"\eta(t, x)", fontsize)),
          linewidth = 2)

    savefig(joinpath(OUT, "bathymetry.pdf"))
end

# Plot of diserpersion relations
function fig_2()
    linewidth = 2
    markersize = 5

    h0 = 1.0
    g = 1.0
    c0 = sqrt(g * h0)

    k = 0.01:0.5:(8 * pi)
    k_zoom = 0.01:0.3:pi
    ylim = (0.0, 1.1)

    omega_euler(k) = sqrt(g * k) * sqrt(tanh(h0 * k))
    c_euler(k) = omega_euler(k) / k
    plot(k, c_euler.(k) ./ c0, label = "Euler", ylim = ylim, xguide = L"k",
         yguide = L"c/c_0", linewidth = linewidth, markershape = :circle,
         markersize = markersize)
    plot!(k_zoom, c_euler.(k_zoom) ./ c0, ylim = (0.54, 1.0),
          inset = bbox(0.35, 0.1, 0.35, 0.3), subplot = 2, legend = nothing,
          linewidth = linewidth, markershape = :circle, markersize = markersize,
          framestyle = :box)

    function plot_dispersion_relation(omega, label, markershape)
        c(k) = omega(k) / k
        plot!(k, c.(k) ./ c0, label = label, linewidth = linewidth,
              markershape = markershape, markersize = markersize)
        plot!(k_zoom, c.(k_zoom) ./ c0, subplot = 2, linewidth = linewidth,
              markershape = markershape, markersize = markersize)
    end

    omega_bbmbbm_(k, d0) = sqrt(g * h0) * k / (1 + 1 / 6 * (d0 * k)^2)
    omega_bbmbbm(k) = omega_bbmbbm_(k, h0)
    plot_dispersion_relation(omega_bbmbbm, "BBM-BBM", :cross)

    alpha_set1 = -1 / 3 * c0 * h0^2
    beta_set1 = 0.0 * h0^3
    gamma_set1 = 0.0 * c0 * h0^3

    alpha_set2 = 0.0004040404040404049 * c0 * h0^2
    beta_set2 = 0.49292929292929294 * h0^3
    gamma_set2 = 0.15707070707070708 * c0 * h0^3

    alpha_set3 = 0.0 * c0 * h0^2
    beta_set3 = 0.27946992481203003 * h0^3
    gamma_set3 = 0.0521077694235589 * c0 * h0^3

    alpha_set4 = 0.0 * c0 * h0^2
    beta_set4 = 0.2308939393939394 * h0^3
    gamma_set4 = 0.04034343434343434 * c0 * h0^3

    function char_equation(alpha, beta, gamma, k)
        a = (1 + beta / h0 * k^2)
        b = (-alpha - beta * alpha / h0 * k^2 - gamma / h0) * k^3
        c = -g * h0 * k^2 + gamma * alpha / h0 * k^6
        omega1 = (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
        #         omega2 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
        return omega1
    end

    omega_set1(k) = char_equation(alpha_set1, beta_set1, gamma_set1, k)
    plot_dispersion_relation(omega_set1, "S.-K. set 1", :rtriangle)

    omega_set2(k) = char_equation(alpha_set2, beta_set2, gamma_set2, k)
    plot_dispersion_relation(omega_set2, "S.-K. set 2", :star5)

    omega_set3(k) = char_equation(alpha_set3, beta_set3, gamma_set3, k)
    plot_dispersion_relation(omega_set3, "S.-K. set 3", :star8)

    omega_set4(k) = char_equation(alpha_set4, beta_set4, gamma_set4, k)
    plot_dispersion_relation(omega_set4, "S.-K. set 4", :diamond)

    # Plot box
    plot!([0.0, pi], [0.54, 0.54], color = :black, label = :none)
    plot!([0.0, pi], [1.0, 1.0], color = :black, label = :none)
    plot!([0.0, 0.0], [0.54, 1.0], color = :black, label = :none)
    plot!([pi, pi], [0.54, 1.0], color = :black, label = :none)

    # Plot connecting lines
    plot!([pi, 6.8], [0.54, 0.629], color = :black, label = :none)
    plot!([pi, 6.8], [1, 1.01], color = :black, label = :none)

    savefig(joinpath(OUT, "dispersion_relations.pdf"))
end

# Chapter 4
# Chaper 4.1 Soliton

const OUT_SOLITON = joinpath(OUT, "soliton")
ispath(OUT_SOLITON) || mkpath(OUT_SOLITON)

# Plot convergence orders for baseline and relaxation
function fig_3()
    tspan = (0.0, 10.0)
    accuracy_orders = [2, 4, 6, 8]
    iters = [4, 4, 4, 3]
    initial_Ns = [128, 128, 128, 128]

    all_Ns = minimum(initial_Ns) * 2 .^ (0:(maximum(iters) - 1))

    linewidth = 2
    markersize = 5
    markershapes = [:circle, :star5, :star8, :rtriangle]
    plot(label = :none, xscale = :log2, yscale = :log10, xlabel = "N", ylim = (1e-5, 1e2),
         ylabel = L"\Vert\eta - \eta_{ana}\Vert_2 + \Vert v - v_{ana}\Vert_2",
         legend = :bottomleft, layout = 2)

    # left subplot: baseline
    for i in 1:length(accuracy_orders)
        Ns = initial_Ns[i] * 2 .^ (0:(iters[i] - 1))
        _, errormatrix = convergence_test("examples/bbm_bbm_1d/bbm_bbm_1d_basic.jl",
                                          iters[i]; N = initial_Ns[i], tspan = tspan,
                                          accuracy_order = accuracy_orders[i])
        # Use sum over all L^2-errors for all variables, i.e. ||η - η_ana||_2 + ||v - v_ana||_2
        l2_err = sum(errormatrix[:l2], dims = 2)
        eocs = log.(l2_err[2:end] ./ l2_err[1:(end - 1)]) ./ log(0.5)
        eoc_mean = round(sum(eocs) / length(eocs), digits = 2)
        plot!(Ns, l2_err, label = "p = $accuracy_order, EOC: $eoc_mean",
              markershape = markershapes[i], linewidth = linewidth, markersize = markersize,
              subplot = 1)
    end
    xticks!(all_Ns, string.(all_Ns), subplot = 1)

    # right subplot: relaxation
    for i in 1:length(accuracy_orders)
        Ns = initial_Ns[i] * 2 .^ (0:(iters[i] - 1))
        _, errormatrix = convergence_test("examples/bbm_bbm_1d/bbm_bbm_1d_relaxation.jl",
                                          iters[i]; N = initial_Ns[i], tspan = tspan,
                                          accuracy_order = accuracy_orders[i])
        # Use sum over all L^2-errors for all variables, i.e. ||η - η_ana||_2 + ||v - v_ana||_2
        l2_err = sum(errormatrix[:l2], dims = 2)
        eocs = log.(l2_err[2:end] ./ l2_err[1:(end - 1)]) ./ log(0.5)
        eoc_mean = round(sum(eocs) / length(eocs), digits = 2)
        plot!(Ns, l2_err, label = "p = $accuracy_order, EOC: $eoc_mean",
              markershape = markershapes[i], linewidth = linewidth, markersize = markersize,
              subplot = 2)
    end
    xticks!(all_Ns, string.(all_Ns), subplot = 2)
    savefig(joinpath(OUT_SOLITON, "orders.pdf"))
end

# Plot errors, change of invariants, and solution at final time for baseline and relaxation
function fig_4_5_6()
    linewidth = 2
    linestyles = [:dash, :dot]

    g = 9.81
    D = 2.0
    c = 5 / 2 * sqrt(g * D)
    x_min = -35.0
    x_max = 35.0
    tspan = (0.0, 50 * (x_max - x_min) / c)
    N = 512
    accuracy_order = 8

    # baseline
    run_example(joinpath(EXAMPLES_DIR_BBMBBM, "bbm_bbm_1d_basic.jl"),
                gravity_constant = g, D = D, coordinates_min = x_min,
                coordinates_max = x_max, tspan = tspan, N = N,
                accuracy_order = accuracy_order)
    p1 = plot(analysis_callback, title = "", label_extension = "baseline",
              linestyles = [:solid :dash :dot],
              linewidth = linewidth, layout = 2, subplot = 1)
    p2 = plot(analysis_callback, title = "", what = (:errors,),
              label_extension = "baseline", linestyle = linestyles[1],
              linewidth = linewidth,
              ylabel = L"\Vert\eta - \eta_{ana}\Vert_2 + \Vert v - v_{ana}\Vert_2",
              exclude = [:conservation_error])
    p3 = plot(semi => sol, label = "baseline", plot_initial = true, plot_bathymetry = false,
              linestyle = linestyles[1], linewidth = linewidth, plot_title = "", title = "",
              ylims = [(-8, 3) (-1, 40)])
    x = DispersiveShallowWater.grid(semi)
    q = DispersiveShallowWater.wrap_array(sol.u[end], semi)
    plot!(p3, x, view(q, 1, :), inset = (1, bbox(0.11, 0.6, 0.35, 0.32)), subplot = 3,
          xlim = (-20, -10),
          ylim = (-0.05, 0.05), legend = nothing, linewidth = linewidth,
          linestyle = linestyles[1],
          color = 2,
          tickfontsize = 5, yticks = [0.04, 0.0, -0.04], xticks = [-20, -15, -10],
          plot_initial = true, plot_bathymetry = false, framestyle = :box)
    q_exact = DispersiveShallowWater.wrap_array(DispersiveShallowWater.compute_coefficients(initial_condition,
                                                                                            tspan[2],
                                                                                            semi),
                                                semi)
    plot!(p3, x, view(q_exact, 1, :), subplot = 3, legend = nothing, linewidth = linewidth,
          linestyle = :solid, color = 1)

    # relaxation
    run_example(joinpath(EXAMPLES_DIR_BBMBBM, "bbm_bbm_1d_relaxation.jl"),
                gravity_constant = g, D = D, coordinates_min = x_min,
                coordinates_max = x_max, tspan = tspan, N = N,
                accuracy_order = accuracy_order)
    plot!(p1, analysis_callback, title = "", label_extension = "relaxation",
          style = [:solid :dash :dot],
          linewidth = linewidth, subplot = 2)
    plot!(p2, analysis_callback, title = "", what = (:errors,),
          label_extension = "relaxation", linestyle = linestyles[2], linewidth = linewidth,
          ylabel = L"\Vert\eta - \eta_{ana}\Vert_2 + \Vert v - v_{ana}\Vert_2",
          exclude = [:conservation_error])
    plot!(p3, semi => sol, plot_bathymetry = false, label = "relaxation",
          linestyle = linestyles[2],
          linewidth = linewidth, plot_title = "", title = "", color = 3)
    x = DispersiveShallowWater.grid(semi)
    q = DispersiveShallowWater.wrap_array(sol.u[end], semi)
    plot!(p3, x, view(q, 1, :), subplot = 3, legend = nothing, linewidth = linewidth,
          linestyle = linestyles[2], color = 3)

    # Plot box
    plot!(p3, [-20, -10], [-0.1, -0.1], color = :black, label = :none)
    plot!(p3, [-20, -10], [0.1, 0.1], color = :black, label = :none)
    plot!(p3, [-20, -20], [-0.1, 0.1], color = :black, label = :none)
    plot!(p3, [-10, -10], [-0.1, 0.1], color = :black, label = :none)

    # Plot connecting lines
    plot!(p3, [-20, -29], [-0.1, -3.6], color = :black, label = :none)
    plot!(p3, [-10, -3.15], [-0.1, -3.6], color = :black, label = :none)

    savefig(p1, joinpath(OUT_SOLITON, "invariants.pdf"))
    savefig(p2, joinpath(OUT_SOLITON, "errors.pdf"))
    savefig(p3, joinpath(OUT_SOLITON, "solution.pdf"))
end

# Chapter 4.3 Dingemans
const OUT_DINGEMANS = joinpath(OUT, "dingemans")
ispath(OUT_DINGEMANS) || mkpath(OUT_DINGEMANS)

# Plot of total waterheight for different models at different points in time
function fig_7()
    linewidth = 3
    fontsize = 16
    linestyles = [:solid, :dash, :dot]

    N = 512
    steps = [100, 200, 300, 500]
    xlims_zoom = [(-25, 0), (5, 30), (20, 45), (-100, -75)]
    ylim_zoom = (0.75, 0.85)
    run_example(joinpath(EXAMPLES_DIR_BBMBBM_VARIABLE,
                         "bbm_bbm_variable_bathymetry_1d_dingemans.jl");
                N = N)
    plot(layout = (2, 2), ylim = (-0.05, 0.86), size = (1200, 800),
         titlefontsize = fontsize)
    for (i, step) in enumerate(steps)
        plot!(semi => sol, step = step, conversion = waterheight_total, label = "BBM-BBM",
              subplot = i, plot_title = "", linewidth = linewidth, legend = :none,
              guidefontsize = fontsize, tickfontsize = fontsize, linestyle = linestyles[1])
        plot!(semi => sol, step = step, inset = (i, bbox(0.1, 0.2, 0.6, 0.5)),
              conversion = waterheight_total, linewidth = linewidth, legend = :none,
              framestyle = :box, xlim = xlims_zoom[i], ylim = ylim_zoom,
              subplot = length(steps) + i, plot_title = "", title = "", xguide = "",
              yguide = "", linestyle = linestyles[1])
    end

    run_example(joinpath(EXAMPLES_DIR_SVAERD_KALISCH, "svaerd_kalisch_1d_dingemans.jl");
                N = N)
    for (i, step) in enumerate(steps)
        plot!(semi => sol, step = step, conversion = waterheight_total,
              label = "Svärd-Kalisch (set 3)", subplot = i, plot_bathymetry = false,
              plot_title = "", linewidth = linewidth, legend = :none,
              guidefontsize = fontsize, tickfontsize = fontsize, color = 2,
              linestyle = linestyles[2])
        plot!(semi => sol, step = step, conversion = waterheight_total,
              linewidth = linewidth, legend = :none, framestyle = :box,
              xlim = xlims_zoom[i], ylim = ylim_zoom, subplot = length(steps) + i,
              plot_title = "", title = "", xguide = "", yguide = "", color = 2,
              linestyle = linestyles[2])
    end

    include("elixir_shallowwater_1d_dingemans.jl")
    for (i, step) in enumerate(steps)
        plot!(PlotData1D(sol.u[step], semi)["H"], label = "Shallow water", subplot = i,
              title = "t = $(round(sol.t[step], digits = 2))", plot_title = "",
              linewidth = linewidth, legend = :none, guidefontsize = fontsize,
              tickfontsize = fontsize, color = 3, linestyle = linestyles[3])
        plot!(PlotData1D(sol.u[step], semi)["H"], linewidth = linewidth, legend = :none,
              framestyle = :box, xlim = xlims_zoom[i], ylim = ylim_zoom,
              subplot = length(steps) + i, plot_title = "", title = "", xguide = "",
              yguide = "", color = 3, linestyle = linestyles[3])
    end

    # dirty hack to have one legend for all subplots
    plot!(subplot = 3, legend_column = 2, bottom_margin = 22 * Plots.mm,
          legend = (0.7, -0.34), legendfontsize = 12)
    plot!(left_margin = 5 * Plots.mm)

    # plot boxes
    for i in 1:length(steps)
        plot!([xlims_zoom[i][1], xlims_zoom[i][2]], [ylim_zoom[1], ylim_zoom[1]],
              color = :black, label = :none, subplot = i, linewidth = 2)
        plot!([xlims_zoom[i][1], xlims_zoom[i][2]], [ylim_zoom[2], ylim_zoom[2]],
              color = :black, label = :none, subplot = i, linewidth = 2)
        plot!([xlims_zoom[i][1], xlims_zoom[i][1]], [ylim_zoom[1], ylim_zoom[2]],
              color = :black, label = :none, subplot = i, linewidth = 2)
        plot!([xlims_zoom[i][2], xlims_zoom[i][2]], [ylim_zoom[1], ylim_zoom[2]],
              color = :black, label = :none, subplot = i, linewidth = 2)
    end
    # plot connecting lines
    upper_corners = [[-119.5, 0.68], [-9.5, 0.68]]
    for i in 1:length(steps)
        plot!([xlims_zoom[i][1], upper_corners[1][1]], [ylim_zoom[1], upper_corners[1][2]],
              color = :black, label = :none, subplot = i, linewidth = 2)
        plot!([xlims_zoom[i][2], upper_corners[2][1]], [ylim_zoom[1], upper_corners[2][2]],
              color = :black, label = :none, subplot = i, linewidth = 2)
    end
    savefig(joinpath(OUT_DINGEMANS, "waterheight_over_time.pdf"))
end

# Plot of total waterheight for Svärd-Kalisch equations at different points in space and different orders of accuracy
function fig_8()
    ylim = (0.75, 0.85)
    yticks = [0.76, 0.78, 0.8, 0.82, 0.84]
    x_values = [3.04, 9.44, 20.04, 26.04, 30.44, 37.04]
    tlims = [
        (15.0, 45.0),
        (19.0, 48.0),
        (25.0, 52.0),
        (30.0, 60.0),
        (33.0, 61.0),
        (35.0, 65.0),
    ]
    plot(layout = (3, 2))

    N = 512
    tspan = (0.0, 70.0)
    saveat = range(tspan..., length = 1000)
    accuracy_orders = [2, 4, 6]
    linestyles = [:solid, :dash, :dot]

    for (i, accuracy_order) in enumerate(accuracy_orders)
        run_example(joinpath(EXAMPLES_DIR_SVAERD_KALISCH, "svaerd_kalisch_1d_dingemans.jl");
                    N = N, tspan = tspan, accuracy_order = accuracy_order, saveat = saveat)
        for (j, x) in enumerate(x_values)
            index = argmin(abs.(DispersiveShallowWater.grid(semi) .- x))
            title = "x = $(round(DispersiveShallowWater.grid(semi)[index], digits = 4))"
            plot!(semi => sol, x, conversion = waterheight_total, subplot = j,
                  xlim = tlims[j], ylim = ylim, plot_title = "", title = title,
                  legend = nothing, yticks = yticks, linewidth = 2, titlefontsize = 10,
                  label = "p = $accuracy_order ", linestyle = linestyles[i])
        end
    end

    plot!(subplot = 5, legend = (0.82, -1.0), legend_column = 3, legendfontsize = 8,
          bottom_margin = 10 * Plots.mm)
    savefig(joinpath(OUT_DINGEMANS, "waterheight_at_x_accuracy_order.pdf"))
end

# Plots of total waterheight for Svärd-Kalisch equations at different points in space and different types of solvers
function fig_9()
    ylim = (0.75, 0.85)
    yticks = [0.76, 0.78, 0.8, 0.82, 0.84]
    x_values = [3.04, 9.44, 20.04, 26.04, 30.44, 37.04]
    tlims = [
        (15.0, 45.0),
        (19.0, 48.0),
        (25.0, 52.0),
        (30.0, 60.0),
        (33.0, 61.0),
        (35.0, 65.0),
    ]
    plot(layout = (3, 2))

    N = 512
    tspan = (0.0, 70.0)
    saveat = range(tspan..., length = 1000)
    accuracy_order = 4
    linestyles = [:solid, :dash, :dot]

    p = 3 # N needs to be divisible by p + 1
    D_legendre = legendre_derivative_operator(-1.0, 1.0, p + 1)
    uniform_mesh = UniformPeriodicMesh1D(coordinates_min, coordinates_max, div(N, p + 1))
    D1 = couple_discontinuously(D_legendre, uniform_mesh)
    D_pl = couple_discontinuously(D_legendre, uniform_mesh, Val(:plus))
    D_min = couple_discontinuously(D_legendre, uniform_mesh, Val(:minus))
    D2 = sparse(D_pl) * sparse(D_min)
    solver_DG = Solver(D1, D2)

    p = 4 # N needs to be divisible by p
    D_legendre = legendre_derivative_operator(-1.0, 1.0, p + 1)
    uniform_mesh = UniformPeriodicMesh1D(coordinates_min, coordinates_max, div(N, p))
    D1 = couple_continuously(D_legendre, uniform_mesh)
    D2_legendre = legendre_second_derivative_operator(-1.0, 1.0, p + 1)
    D2 = couple_continuously(D2_legendre, uniform_mesh)
    solver_CG = Solver(D1, D2)

    solvers = [solver_DG, :none, solver_CG]
    labels = ["DG ", "FD ", "CG "]

    for (i, solver) in enumerate(solvers)
        if solver == :none
            run_example(joinpath(EXAMPLES_DIR_SVAERD_KALISCH,
                                 "svaerd_kalisch_1d_dingemans.jl");
                        N = N, tspan = tspan, accuracy_order = accuracy_order,
                        saveat = saveat)
        else
            run_example(joinpath(EXAMPLES_DIR_SVAERD_KALISCH,
                                 "svaerd_kalisch_1d_dingemans.jl");
                        N = N, tspan = tspan, accuracy_order = accuracy_order,
                        saveat = saveat, solver = solvers[i])
        end
        for (j, x) in enumerate(x_values)
            index = argmin(abs.(DispersiveShallowWater.grid(semi) .- x))
            title = "x = $(round(DispersiveShallowWater.grid(semi)[index], digits = 4))"
            plot!(semi => sol, x, conversion = waterheight_total, subplot = j,
                  xlim = tlims[j], ylim = ylim, plot_title = "", title = title,
                  legend = nothing, yticks = yticks, linewidth = 2, titlefontsize = 10,
                  label = labels[i], linestyle = linestyles[i])
        end
    end

    plot!(subplot = 5, legend = (0.86, -1.0), legend_column = 3, legendfontsize = 8,
          bottom_margin = 10 * Plots.mm)
    savefig(joinpath(OUT_DINGEMANS, "waterheight_at_x_solver_types.pdf"))
end

function fig_10_11()
    ylim = (0.75, 0.85)
    yticks = [0.76, 0.78, 0.8, 0.82, 0.84]
    x_values = [3.04, 9.44, 20.04, 26.04, 30.44, 37.04]
    tlims = [
        (15.0, 45.0),
        (19.0, 48.0),
        (25.0, 52.0),
        (30.0, 60.0),
        (33.0, 61.0),
        (35.0, 65.0),
    ]
    p1 = plot(layout = (3, 2))

    N = 512
    tspan = (0.0, 70.0)
    saveat = range(tspan..., length = 1000)
    accuracy_order = 4
    linestyles = [:solid, :dash, :dot]
    linewidth = 2
    titlefontsize = 10

    labels = ["EC baseline ", "EC relaxation ", "ED "]

    function plot_at_x(semi, sol, i)
        for (j, x) in enumerate(x_values)
            index = argmin(abs.(DispersiveShallowWater.grid(semi) .- x))
            title = "x = $(round(DispersiveShallowWater.grid(semi)[index], digits = 4))"
            plot!(p1, semi => sol, x, conversion = waterheight_total, subplot = j,
                  xlim = tlims[j], ylim = ylim, plot_title = "", title = title,
                  legend = nothing, yticks = yticks, linewidth = linewidth,
                  titlefontsize = titlefontsize,
                  label = labels[i], linestyle = linestyles[i])
        end
    end

    run_example(joinpath(EXAMPLES_DIR_SVAERD_KALISCH, "svaerd_kalisch_1d_dingemans.jl");
                N = N, tspan = tspan, accuracy_order = accuracy_order, saveat = saveat)
    plot_at_x(semi, sol, 1)
    p2 = plot(analysis_callback, title = labels[1], legend = :none,
              linestyles = [:solid :dash :dot],
              linewidth = linewidth, layout = 3, subplot = 1, titlefontsize = titlefontsize)

    run_example(joinpath(EXAMPLES_DIR_SVAERD_KALISCH,
                         "svaerd_kalisch_1d_dingemans_relaxation.jl");
                N = N, tspan = tspan, accuracy_order = accuracy_order, saveat = saveat)
    plot_at_x(semi, sol, 2)
    plot!(p2, analysis_callback, title = labels[2], legend = :none,
          linestyles = [:solid :dash :dot],
          linewidth = linewidth, subplot = 2, titlefontsize = titlefontsize)

    run_example(joinpath(EXAMPLES_DIR_SVAERD_KALISCH,
                         "svaerd_kalisch_1d_dingemans_upwind.jl");
                N = N, tspan = tspan, accuracy_order = accuracy_order, saveat = saveat)
    plot_at_x(semi, sol, 3)
    plot!(p2, analysis_callback, title = labels[3], legend = :none,
          linestyles = [:solid :dash :dot],
          linewidth = linewidth, subplot = 3, titlefontsize = titlefontsize)

    plot!(p1, subplot = 5, legend = (0.86, -1.0), legend_column = 3, legendfontsize = 8,
          bottom_margin = 10 * Plots.mm)
    plot!(p2, subplot = 3, legend = (1.3, 0.6), legendfontsize = 8)
    savefig(p1, joinpath(OUT_DINGEMANS, "waterheight_at_x_ec.pdf"))
    savefig(p2, joinpath(OUT_DINGEMANS, "invariants_ec.pdf"))
end

fig_1()
fig_2()
fig_3()
fig_4_5_6()
fig_7()
fig_8()
fig_9()
fig_10_11()
