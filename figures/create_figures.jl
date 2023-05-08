using DispersiveShallowWater
using Plots

EXAMPLES_DIR = joinpath(pathof(DispersiveShallowWater) |> dirname |> dirname, "examples")

function create_gif(filename)
  include(joinpath(EXAMPLES_DIR, filename))
  x = DispersiveShallowWater.grid(semi)
  anim = @animate for step in 1:length(sol.u)
    t = sol.t[step]
    plot(x, sol.u[step][1, :], legend=true, ylim=(-4.0, 2.0), label="approximation", title="time t=$(round(t, digits=5))")
    analytical_sol = DispersiveShallowWater.compute_coefficients(initial_condition, solver, t, equations, mesh)
    plot!(x, analytical_sol[1, :], legend=true, label="analytical")
  end
  isdir("out") || mkdir("out")
  gif(anim, "out/solution_" * splitext(filename)[1] * ".gif", fps = 50)
end

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using periodic SBP operators
create_gif("bbm_bbm_1d_basic.jl")

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using discontinuously coupled Legendre SBP operators
create_gif("bbm_bbm_1d_dg.jl")
