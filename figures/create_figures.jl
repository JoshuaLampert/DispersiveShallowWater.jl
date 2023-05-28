using DispersiveShallowWater
using Plots

EXAMPLES_DIR = pkgdir(DispersiveShallowWater, "examples")

function create_gif(filename; ylims=:auto)
  include(joinpath(EXAMPLES_DIR, filename))
  x = DispersiveShallowWater.grid(semi)
  anim = @animate for step in 1:length(sol.u)
    t = sol.t[step]
    plot(x, view(sol.u[step], 1, :), legend=true, ylims=ylims, label="approximation", title="time t=$(round(t, digits=5))")
    analytical_sol = DispersiveShallowWater.compute_coefficients(initial_condition, solver, t, equations, mesh)
    plot!(x, view(analytical_sol, 1, :), legend=true, label="analytical")
  end
  isdir("out") || mkdir("out")
  gif(anim, "out/solution_" * splitext(filename)[1] * ".gif", fps = 50)
end

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using periodic SBP operators
create_gif("bbm_bbm_1d_basic.jl", ylims=(-8, 4))

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using discontinuously coupled Legendre SBP operators
create_gif("bbm_bbm_1d_dg.jl", ylims=(-4, 2))
