module TestBBMBBMVariableBathymetry1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "svaerd_kalisch_1d")

@testset "SvaerdKalisch1D" begin
    @trixi_testset "svaerd_kalisch_1d_dingemans" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_dingemans.jl"),
                            tspan=(0.0, 1.0),
                            N=512,
                            l2=[0.22650314162825236 0.7326320951882187 0.0],
                            linf=[0.03627318072294938 0.11690541155100909 0.0],
                            cons_error=[3.979039320256561e-13 7.020733524160366e-5 0.0],
                            change_waterheight=-3.979039320256561e-13,
                            change_momentum=-2.0220130084180887e-9,
                            change_entropy=-0.0011575500814160478,
                            change_entropy_modified=-0.0005787808884178958)
    end

    @trixi_testset "svaerd_kalisch_1d_dingemans_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_dingemans_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            N=512,
                            l2=[0.22646027715939582 0.7325839439787915 0.0],
                            linf=[0.03626762302524622 0.11692754047634027 0.0],
                            cons_error=[3.979039320256561e-13 7.038992525406615e-5 0.0],
                            change_waterheight=-3.979039320256561e-13,
                            change_entropy=-0.0005989344795125362,
                            change_entropy_modified=0.0)
    end

    @trixi_testset "svaerd_kalisch_1d_well_balanced" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_well_balanced.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.0 1.135448143093612e-14 0.0],
                            linf=[0.0 8.133477278069499e-15 0.0],
                            cons_error=[0.0 1.6056589579882354e-14 0.0],
                            change_waterheight=0.0,
                            change_momentum=1.5679986322667355e-14,
                            change_entropy=0.0,
                            lake_at_rest=0.0)
    end
end

end # module
