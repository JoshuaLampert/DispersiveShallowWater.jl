using DispersiveShallowWater
using LaTeXStrings
using Plots

function plot_elixir(filename; ylims_gif = nothing, ylims_x = nothing, x_values = [],
                     tlims = fill(:auto, length(x_values)), kwargs...)
    trixi_include(joinpath(examples_dir(), filename); kwargs...)
    elixirname = splitext(basename(filename))[1]
    outdir = joinpath("out", dirname(filename), elixirname)
    ispath(outdir) || mkpath(outdir)

    # Plot solution
    anim = @animate for step in 1:length(sol.u)
        plot(semi => sol, plot_initial = true, step = step, yli = ylims_gif)
    end
    gif(anim, joinpath(outdir, "solution.gif"), fps = 25)

    # Plot error in invariants
    plot(analysis_callback)
    savefig(joinpath(outdir, "invariants.pdf"))

    # Plot at different x values over time
    @assert size(x_values) == size(tlims)
    for (i, x) in enumerate(x_values)
        plot(semi => sol, x, xlim = tlims[i], yli = ylims_x)
        savefig(joinpath(outdir, "solution_at_x_" * string(x) * ".pdf"))
    end
end

EXAMPLES_DIR_BBMBBM = "bbm_bbm_1d"
EXAMPLES_DIR_BBMBBM_VARIABLE = "bbm_bbm_variable_bathymetry_1d"
EXAMPLES_DIR_SVAERD_KALISCH = "svaerd_kalisch_1d"

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using periodic SBP operators
plot_elixir(joinpath(EXAMPLES_DIR_BBMBBM, "bbm_bbm_1d_basic.jl");
            ylims_gif = [(-8, 4), :auto], tspan = (0.0, 50.0))

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using discontinuously coupled Legendre SBP operators
plot_elixir(joinpath(EXAMPLES_DIR_BBMBBM, "bbm_bbm_1d_dg.jl");
            ylims_gif = [(-4, 2), :auto])

###############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using periodic SBP operators and relaxation, is energy-conservative
plot_elixir(joinpath(EXAMPLES_DIR_BBMBBM, "bbm_bbm_1d_relaxation.jl");
            ylims_gif = [(-8, 4), (-10, 30)],
            tspan = (0.0, 30.0))

##############################################################################
# Travelling wave solution for one-dimensional BBM-BBM equations with periodic boundary conditions
# using periodic SBP operators. Uses the BBM-BBM equations with variable bathymetry, but sets the bathymetry
# as a constant. Should give the same result as "bbm_bbm_1d_basic.jl"
plot_elixir(joinpath(EXAMPLES_DIR_BBMBBM_VARIABLE,
                     "bbm_bbm_variable_bathymetry_1d_basic.jl");
            ylims_gif = [(-8, 4), :auto],
            tspan = (0.0, 50.0))

###############################################################################
# One-dimensional BBM-BBM equations with a Gaussian bump as initial condition for the water height
# and initially still water. The bathymetry is a sine function. Relaxation is used, so the solution
# is energy-conservative. Uses periodic finite difference SBP operators.
plot_elixir(joinpath(EXAMPLES_DIR_BBMBBM_VARIABLE,
                     "bbm_bbm_variable_bathymetry_1d_relaxation.jl");
            ylims_gif = [(-1.5, 6.0), (-10.0, 10.0)],
            tspan = (0.0, 10.0))

###############################################################################
# One-dimensional BBM-BBM equations with a Gaussian bump as initial condition for the water height
# and initially still water. The bathymetry is a sine function. Relaxation is used, so the solution
# is energy-conservative. Uses upwind discontinuously coupled SBP operators.
plot_elixir(joinpath(EXAMPLES_DIR_BBMBBM_VARIABLE,
                     "bbm_bbm_variable_bathymetry_1d_dg_upwind_relaxation.jl");
            ylims_gif = [(-1.5, 6.0), (-10.0, 10.0)],
            tspan = (0.0, 10.0))

###############################################################################
# One-dimensional BBM-BBM equations with a Gaussian bump as initial condition for the water height
# and initially still water. The bathymetry is a sine function. Relaxation is used, so the solution
# is energy-conservative. Uses periodic finite difference discontinuously coupled SBP operators.
plot_elixir(joinpath(EXAMPLES_DIR_BBMBBM_VARIABLE,
                     "bbm_bbm_variable_bathymetry_1d_upwind_relaxation.jl");
            ylims_gif = [(-1.5, 6.0), (-10.0, 10.0)],
            tspan = (0.0, 10.0))

###############################################################################
# One-dimensional BBM-BBM equations with a constant water height
# and initially still water. The bathymetry is discontinuous. Relaxation is used, so the solution
# is energy-conservative. Uses periodic finite difference SBP operators. The solution should be
# (exactly) constant in time.
plot_elixir(joinpath(EXAMPLES_DIR_BBMBBM_VARIABLE,
                     "bbm_bbm_variable_bathymetry_1d_well_balanced.jl");
            ylims_gif = [(2.0 - 1e-3, 2.0 + 1e-3), (-1e-3, 1e-3)],
            tspan = (0.0, 10.0))

##############################################################################
# One-dimensional BBM-BBM equations with initial condition that models
# a wave make. This setup comes from experiments by W. M. Dingemans.
plot_elixir(joinpath(EXAMPLES_DIR_BBMBBM_VARIABLE,
                     "bbm_bbm_variable_bathymetry_1d_dingemans.jl");
            ylims_gif = [(-0.1, 0.9), (-0.3, 0.3)],
            ylims_x = [:auto, :auto],
            x_values = [3.04, 9.44, 20.04, 26.04, 30.44, 37.04],
            tlims = [
                (15.0, 45.0),
                (19.0, 48.0),
                (25.0, 52.0),
                (30.0, 60.0),
                (33.0, 61.0),
                (35.0, 65.0),
            ],
            tspan = (0.0, 70.0))

###############################################################################
# One-dimensional equations from Svärd and Kalisch with initial condition that models
# a wave make. This setup comes from experiments by W. M. Dingemans.
plot_elixir(joinpath(EXAMPLES_DIR_SVAERD_KALISCH,
                     "svaerd_kalisch_1d_dingemans.jl");
            ylims_gif = [(-0.1, 0.9), (-0.3, 0.3)],
            ylims_x = [:auto, :auto],
            x_values = [3.04, 9.44, 20.04, 26.04, 30.44, 37.04],
            tlims = [
                (15.0, 45.0),
                (19.0, 48.0),
                (25.0, 52.0),
                (30.0, 60.0),
                (33.0, 61.0),
                (35.0, 65.0),
            ],
            tspan = (0.0, 70.0))

###############################################################################
# One-dimensional Svärd-Kalisch equations with a constant water height
# and initially still water. The bathymetry is discontinuous. Relaxation is used, so the solution
# is energy-conservative. Uses periodic finite difference SBP operators. The solution should be
# (exactly) constant in time.
plot_elixir(joinpath(EXAMPLES_DIR_SVAERD_KALISCH,
                     "svaerd_kalisch_1d_well_balanced.jl");
            ylims = [(2.0 - 1e-3, 2.0 + 1e-3), (-1e-3, 1e-3)],
            tspan = (0.0, 10.0))
