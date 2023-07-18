module TestShallowWater1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "shallow_water_1d")

@testset "ShallowWater1D" begin
    @trixi_testset "shallow_water_1d_basic" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "shallow_water_1d_basic.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.22941778911937488 0.3489735105452126 0.0],
                            linf=[0.3029397223767515 0.5993746621524781 0.0],
                            cons_error=[1.7763568394002505e-15 8.315830662963819e-16 0.0],
                            change_waterheight=-1.7763568394002505e-15,
                            change_momentum=5.640019701269594e-16,
                            change_entropy=-3.819764415879945e-7)
    end

    @trixi_testset "shallow_water_1d_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "shallow_water_1d_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.2294178346554772 0.34897380752773605 0.0],
                            linf=[0.3029407021084475 0.5993746890244087 0.0],
                            cons_error=[2.6645352591003757e-15 1.145568015448184e-15 0.0],
                            change_waterheight=-2.6645352591003757e-15,
                            change_momentum=1.3960187172923355e-15,
                            change_entropy=6.394884621840902e-14)
    end

    @trixi_testset "shallow_water_1d_well_balanced" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "shallow_water_1d_well_balanced.jl"),
                            tspan=(0.0, 1.0),
                            l2=[1.2363565604017289e-14 4.0940746032364917e-14 0.0],
                            linf=[2.731148640577885e-14 8.315557448991739e-14 0.0],
                            cons_error=[1.7763568394002505e-15 2.5630855292351704e-14 0.0],
                            change_waterheight=-1.7763568394002505e-15,
                            change_momentum=2.4107318306126452e-14,
                            change_entropy=-2.7000623958883807e-13,
                            lake_at_rest=1.398318324894359e-14)
    end
end

end # module
