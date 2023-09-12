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
                            l2=[0.22647683207160013 0.7325945009751711 0.0],
                            linf=[0.03634466451241691 0.11720005557202748 0.0],
                            cons_error=[3.979039320256561e-13 0.00010706168405306206 0.0],
                            change_waterheight=-3.979039320256561e-13,
                            change_momentum=-2.0220130084180887e-9,
                            change_entropy=-0.001319708553410237,
                            change_entropy_modified=-0.0007363973352312314)
    end

    @trixi_testset "svaerd_kalisch_1d_dingemans_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_dingemans_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            N=512,
                            l2=[0.22648917549343248 0.7327156628336414 0.0],
                            linf=[0.03634866178910967 0.117250088972775337 0.0],
                            cons_error=[3.979039320256561e-13 0.00010728668392322323 0.0],
                            change_waterheight=-3.979039320256561e-13,
                            change_entropy=-0.0006080247380850778,
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
