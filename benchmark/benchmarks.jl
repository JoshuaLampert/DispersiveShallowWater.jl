using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using BenchmarkTools
using DispersiveShallowWater

const SUITE = BenchmarkGroup()

elixirs = [joinpath(examples_dir(), "bbm_1d", "bbm_1d_basic.jl"),
    joinpath(examples_dir(), "bbm_1d", "bbm_1d_fourier.jl"),
    joinpath(examples_dir(), "bbm_bbm_1d", "bbm_bbm_1d_dg.jl"),
    joinpath(examples_dir(), "bbm_bbm_1d", "bbm_bbm_1d_relaxation.jl"),
    joinpath(examples_dir(), "bbm_bbm_1d", "bbm_bbm_1d_upwind_relaxation.jl"),
    joinpath(examples_dir(), "bbm_bbm_1d", "bbm_bbm_1d_basic_reflecting.jl"),
    joinpath(examples_dir(), "hyperbolic_serre_green_naghdi_1d",
             "hyperbolic_serre_green_naghdi_dingemans.jl"),
    joinpath(examples_dir(), "serre_green_naghdi_1d",
             "serre_green_naghdi_well_balanced.jl"),
    joinpath(examples_dir(), "svaerd_kalisch_1d",
             "svaerd_kalisch_1d_dingemans_relaxation.jl")]

for elixir in elixirs
    benchname = joinpath(basename(dirname(elixir)), basename(elixir))
    println("Running $benchname...")
    # use include instead of trixi_include(elixir, tspan = (0.0, 1e-10)) because AirspeedVelocity.jl has namespace issues
    redirect_stdout(devnull) do
        include(elixir)
    end
    SUITE[benchname] = @benchmarkable DispersiveShallowWater.rhs!($(similar(sol.u[end])),
                                                                  $(copy(sol.u[end])),
                                                                  $(semi), $(first(tspan)))
end
