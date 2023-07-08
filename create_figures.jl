using DispersiveShallowWater
using LaTeXStrings
using Plots

function plot_gif_invariants(filename; ylims = nothing, kwargs...)
    trixi_include(joinpath(examples_dir(), filename); kwargs...)
    elixirname = splitext(basename(filename))[1]
    outdir = joinpath("out", dirname(filename))
    ispath(outdir) || mkpath(outdir)

    # Plot solution
    anim = @animate for step in 1:length(sol.u)
        plot(semi => sol, plot_initial = true, step = step, yli = ylims)
    end
    gif(anim, joinpath(outdir, "solution_" * elixirname * ".gif"), fps = 25)

    # Plot error in invariants
    plot(analysis_callback)
    savefig(joinpath(outdir, "invariants_" * elixirname * ".pdf"))
end

EXAMPLES_DIR_BBMBBM = "bbm_bbm_1d"
EXAMPLES_DIR_BBMBBM_VARIABLE = "bbm_bbm_variable_bathymetry_1d"

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using periodic SBP operators
plot_gif_invariants(joinpath(EXAMPLES_DIR_BBMBBM, "bbm_bbm_1d_basic.jl");
                    ylims = [(-8, 4), :auto], tspan = (0.0, 50.0))

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using discontinuously coupled Legendre SBP operators
plot_gif_invariants(joinpath(EXAMPLES_DIR_BBMBBM, "bbm_bbm_1d_dg.jl"); ylims = [(-4, 2), :auto])

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using periodic SBP operators and relaxation, is energy-conservative
plot_gif_invariants(joinpath(EXAMPLES_DIR_BBMBBM, "bbm_bbm_1d_relaxation.jl");
                    ylims = [(-8, 4), (-10, 30)],
                    tspan = (0.0, 30.0))

##############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using periodic SBP operators. Uses the BBM-BBM equations with variable bathymetry, but sets the bathymetry
# as a constant. Should give the same result as "bbm_bbm_1d_basic.jl"
plot_gif_invariants(joinpath(EXAMPLES_DIR_BBMBBM_VARIABLE,
                             "bbm_bbm_variable_bathymetry_1d_basic.jl");
                    ylims = [(-8, 4), :auto],
                    tspan = (0.0, 50.0))

###############################################################################
# One-dimensional BBM-BBM equations with a Gaussian bump as initial condition for the water height
# and initially still water. The bathymetry is a sine function. Relaxation is used, so the solution
# is energy-conservative. Uses periodic finite difference SBP operators.
plot_gif_invariants(joinpath(EXAMPLES_DIR_BBMBBM_VARIABLE,
                             "bbm_bbm_variable_bathymetry_1d_relaxation.jl");
                    ylims = [(-1.5, 6.0), (-10.0, 10.0)],
                    tspan = (0.0, 10.0))

###############################################################################
# One-dimensional BBM-BBM equations with a Gaussian bump as initial condition for the water height
# and initially still water. The bathymetry is a sine function. Relaxation is used, so the solution
# is energy-conservative. Uses upwind discontinuously coupled SBP operators.
plot_gif_invariants(joinpath(EXAMPLES_DIR_BBMBBM_VARIABLE,
                             "bbm_bbm_variable_bathymetry_1d_dg_upwind_relaxation.jl");
                    ylims = [(-1.5, 6.0), (-10.0, 10.0)],
                    tspan = (0.0, 10.0))

###############################################################################
# One-dimensional BBM-BBM equations with a Gaussian bump as initial condition for the water height
# and initially still water. The bathymetry is a sine function. Relaxation is used, so the solution
# is energy-conservative. Uses periodic finite difference discontinuously coupled SBP operators.
plot_gif_invariants(joinpath(EXAMPLES_DIR_BBMBBM_VARIABLE,
                             "bbm_bbm_variable_bathymetry_1d_upwind_relaxation.jl");
                    ylims = [(-1.5, 6.0), (-10.0, 10.0)],
                    tspan = (0.0, 10.0))

###############################################################################
# One-dimensional BBM-BBM equations with a constant water height
# and initially still water. The bathymetry is discontinuous. Relaxation is used, so the solution
# is energy-conservative. Uses periodic finite difference SBP operators. The solution should be
# (exactly) constant in time.
plot_gif_invariants(joinpath(EXAMPLES_DIR_BBMBBM_VARIABLE,
                             "bbm_bbm_variable_bathymetry_1d_well_balanced.jl");
                    ylims = [(2.0 - 1e-3, 2.0 + 1e-3), (-1e-3, 1e-3)],
                    tspan = (0.0, 10.0))
