module TestBBMBBMVariableBathymetry1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "svaerd_kalisch_1d")

@testset "SvaerdKalisch1D" begin
    @trixi_testset "svaerd_kalisch_1d_basic" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_basic.jl"),
                            tspan=(0.0, 1.0),
                            l2=[1.0 1.0 1.0],
                            linf=[1.0 1.0 1.0],
                            cons_error=[1.0 1.0 1.0],
                            change_waterheight=1.0,
                            change_momentum=1.0,
                            change_entropy=1.0)
    end

    @trixi_testset "svaerd_kalisch_1d_dingemans" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_dingemans.jl"),
                            tspan=(0.0, 1.0),
                            l2=[1.0 1.0 1.0],
                            linf=[1.0 1.0 1.0],
                            cons_error=[1.0 1.0 1.0],
                            change_waterheight=1.0,
                            change_momentum=1.0,
                            change_entropy=1.0)
    end

    @trixi_testset "svaerd_kalisch_1d_well_balanced" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_well_balanced.jl"),
                            tspan=(0.0, 1.0),
                            l2=[1.0 1.0 1.0],
                            linf=[1.0 1.0 1.0],
                            cons_error=[1.0 1.0 1.0],
                            change_waterheight=1.0,
                            change_momentum=1.0,
                            change_entropy=1.0,
                            lake_at_rest_error=1.0)
    end
end

end # module
