using DispersiveShallowWater
using LaTeXStrings
using Plots

function plot_gif_invariants(filename; ylims_eta = :auto, ylims_v = :auto, kwargs...)
  trixi_include(joinpath(examples_dir(), filename); kwargs...)
  x = DispersiveShallowWater.grid(semi)

  # Plot solution for eta
  anim = @animate for step in 1:length(sol.u)
    u = DispersiveShallowWater.wrap_array(sol.u[step], semi)
    t = sol.t[step]
    plot(x, view(u, 1, :), legend = true, ylims = ylims_eta, label = "approximation η",
         xlabel = L"x", ylabel = L"\eta", title = "time t=$(round(t, digits=5))")
    u_ode_exact = DispersiveShallowWater.compute_coefficients(initial_condition, t, semi)
    u_exact = DispersiveShallowWater.wrap_array(u_ode_exact, semi)
    plot!(x, view(u_exact, 1, :), legend = true, label = "analytical η")
  end
  isdir("out") || mkdir("out")
  gif(anim, "out/solution_eta_" * splitext(filename)[1] * ".gif", fps = 25)

  # Plot solution for v
  anim = @animate for step in 1:length(sol.u)
    t = sol.t[step]
    u = DispersiveShallowWater.wrap_array(sol.u[step], semi)
    plot(x, view(u, 2, :), legend = true, ylims = ylims_v, label = "approximation v",
         xlabel = L"x", ylabel = L"v", title = "time t=$(round(t, digits=5))")
    u_ode_exact = DispersiveShallowWater.compute_coefficients(initial_condition, t, semi)
    u_exact = DispersiveShallowWater.wrap_array(u_ode_exact, semi)
    plot!(x, view(u_exact, 2, :), legend = true, label = "analytical v")
  end
  gif(anim, "out/solution_v_" * splitext(filename)[1] * ".gif", fps = 25)

  # Plot error in invariants
  tstops = DispersiveShallowWater.tstops(analysis_callback)
  integrals = DispersiveShallowWater.integrals(analysis_callback)
  labels = integral_names(analysis_callback)
  plot(xlabel = "t", ylabel = "change in invariant")
  for i in 1:size(integrals, 1)
    plot!(tstops, view(integrals, i, :) .- view(integrals, i, 1), label = string(labels[i]))
  end
  savefig("out/invariants_" * splitext(filename)[1] * ".pdf")
end

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using periodic SBP operators
plot_gif_invariants("bbm_bbm_1d_basic.jl"; ylims_eta = (-8, 4), tspan = (0.0, 50.0))

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using discontinuously coupled Legendre SBP operators
plot_gif_invariants("bbm_bbm_1d_dg.jl"; ylims_eta = (-4, 2))

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using periodic SBP operators and relaxation, is energy-conservative
plot_gif_invariants("bbm_bbm_1d_relaxation.jl"; ylims_eta = (-8, 4), ylims_v = (-10, 30),
                    tspan = (0.0, 30.0))

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using periodic SBP operators. Uses the BBM-BBM equations with variable bathymetry, but sets the bathymetry
# as a constant. Should give the same result as "bbm_bbm_1d_basic.jl"
plot_gif_invariants("bbm_bbm_variable_bathymetry_1d_basic.jl";
                    ylims_eta = (-8, 4),
                    tspan = (0.0, 50.0))

###############################################################################
# One-dimensional BBM-BBM equations with a Gaussian bump as initial condition for the water height
# and initially still water. The bathymetry is a sine function. Relaxation is used, so the solution
# is energy-conservative. Uses periodic finite difference SBP operators.
plot_gif_invariants("bbm_bbm_variable_bathymetry_1d_relaxation.jl"; ylims_eta = (0, 6),
                    ylims_v = (-10.0, 10.0),
                    tspan = (0.0, 10.0))

###############################################################################
# One-dimensional BBM-BBM equations with a Gaussian bump as initial condition for the water height
# and initially still water. The bathymetry is a sine function. Relaxation is used, so the solution
# is energy-conservative. Uses upwind discontinuously coupled SBP operators.
plot_gif_invariants("bbm_bbm_variable_bathymetry_1d_dg_relaxation.jl"; ylims_eta = (0, 6),
                    ylims_v = (-10.0, 10.0),
                    tspan = (0.0, 10.0))

###############################################################################
# One-dimensional BBM-BBM equations with a constant water height
# and initially still water. The bathymetry is discontinuous. Relaxation is used, so the solution
# is energy-conservative. Uses periodic finite difference SBP operators. The solution should be
# (exactly) constant in time.
plot_gif_invariants("bbm_bbm_variable_bathymetry_1d_well_balanced.jl";
                    ylims_eta = (2.0 - 1e-3, 2.0 + 1e-3),
                    ylims_v = (-1e-3, 1e-3),
                    tspan = (0.0, 10.0))
