using DispersiveShallowWater
using Plots
using LaTeXStrings

const OUT = "out/"
ispath(OUT) || mkpath(OUT)
const EXAMPLES_DIR_BBMBBM = "bbm_bbm_1d"

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
          linewidth = linewidth, markershape = :circle, markersize = markersize, framestyle = :box)

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

    g = 9.81
    D = 2.0
    c = 5 / 2 * sqrt(g * D)
    x_min = -35.0
    x_max = 35.0
    tspan = (0.0, 50 * (x_max - x_min) / c)
    N = 512
    accuracy_order = 8

    # baseline
    trixi_include(joinpath(examples_dir(), EXAMPLES_DIR_BBMBBM, "bbm_bbm_1d_basic.jl"),
                  gravity_constant = g, D = D, coordinates_min = x_min,
                  coordinates_max = x_max, tspan = tspan, N = N,
                  accuracy_order = accuracy_order)
    p1 = plot(analysis_callback, title = "", label_extension = "baseline", style = :auto,
              linewidth = linewidth, layout = 2, subplot = 1)
    p2 = plot(analysis_callback, title = "", what = (:errors,),
              label_extension = "baseline", linestyle = :dash, linewidth = linewidth,
              ylabel = L"\Vert\eta - \eta_{ana}\Vert_2 + \Vert v - v_{ana}\Vert_2",
              exclude = [:conservation_error])
    p3 = plot(semi => sol, label = "baseline", plot_initial = true, linestyle = :dash,
              linewidth = linewidth, plot_title = "", title = "", ylims = [(-8, 3) (-1, 40)])
    x = DispersiveShallowWater.grid(semi)
    q = DispersiveShallowWater.wrap_array(sol.u[end], semi)
    plot!(p3, x, view(q, 1, :), inset = (1, bbox(0.11, 0.6, 0.35, 0.32)), subplot = 3, xlim = (-20, -10),
          ylim = (-0.05, 0.05), legend = nothing, linewidth = linewidth, linestyle = :dash, color = 2,
          tickfontsize = 5, yticks = [0.04, 0.0, -0.04], xticks = [-20, -15, -10], framestyle = :box)
    q_exact = DispersiveShallowWater.wrap_array(DispersiveShallowWater.compute_coefficients(initial_condition, tspan[2], semi), semi)
    plot!(p3, x, view(q_exact, 1, :), subplot = 3, legend = nothing, linewidth = linewidth,
          linestyle = :dot, color = 1)
    # Plot box
    plot!(p3, [-20, -10], [-0.1, -0.1], color = :black, label = :none)
    plot!(p3, [-20, -10], [0.1, 0.1], color = :black, label = :none)
    plot!(p3, [-20, -20], [-0.1, 0.1], color = :black, label = :none)
    plot!(p3, [-10, -10], [-0.1, 0.1], color = :black, label = :none)

    # Plot connecting lines
    plot!(p3, [-20, -29], [-0.1, -3.6], color = :black, label = :none)
    plot!(p3, [-10, -3.15], [-0.1, -3.6], color = :black, label = :none)

    # relaxation
    trixi_include(joinpath(examples_dir(), EXAMPLES_DIR_BBMBBM, "bbm_bbm_1d_relaxation.jl"),
                  gravity_constant = g, D = D, coordinates_min = x_min,
                  coordinates_max = x_max, tspan = tspan, N = N,
                  accuracy_order = accuracy_order)
    plot!(p1, analysis_callback, title = "", label_extension = "relaxation", style = :auto,
          linewidth = linewidth, subplot = 2)
    plot!(p2, analysis_callback, title = "", what = (:errors,),
          label_extension = "relaxation", linestyle = :dot, linewidth = linewidth,
          ylabel = L"\Vert\eta - \eta_{ana}\Vert_2 + \Vert v - v_{ana}\Vert_2",
          exclude = [:conservation_error])
    plot!(p3, semi => sol, plot_bathymetry = false, label = "relaxation", linestyle = :dot,
          linewidth = linewidth, plot_title = "", title = "", color = 3)
    x = DispersiveShallowWater.grid(semi)
    q = DispersiveShallowWater.wrap_array(sol.u[end], semi)
    plot!(p3, x, view(q, 1, :), subplot = 3, legend = nothing, linewidth = linewidth,
          linestyle = :dot, color = 3)

    savefig(p1, joinpath(OUT_SOLITON, "invariants.pdf"))
    savefig(p2, joinpath(OUT_SOLITON, "errors.pdf"))
    savefig(p3, joinpath(OUT_SOLITON, "solution.pdf"))
end

fig_1()
fig_2()
fig_3()
fig_4_5_6()
