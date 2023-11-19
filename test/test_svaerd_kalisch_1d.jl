module TestBBMBBMVariableBathymetry1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "svaerd_kalisch_1d")

@testset "SvaerdKalisch1D" begin
    @trixi_testset "svaerd_kalisch_1d_manufactured" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_manufactured.jl"),
                            tspan=(0.0, 0.1),
                            l2=[7.467887060263923e-5 2.7796353838948894e-8 0.0],
                            linf=[0.0001613144395267163 4.344495230235168e-8 0.0],
                            cons_error=[2.3635607360183997e-16 8.084235123776567e-10 0.0],
                            change_waterheight=-2.3635607360183997e-16,
                            change_entropy=0.1342289500320556)
    end

    @trixi_testset "svaerd_kalisch_1d_dingemans" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_dingemans.jl"),
                            tspan=(0.0, 1.0),
                            N=512,
                            l2=[0.22796106962338855 0.7519327063662515 0.0],
                            linf=[0.036708347831218346 0.12141172207472928 0.0],
                            cons_error=[3.979039320256561e-13 4.937137540373564e-5 0.0],
                            change_waterheight=-3.979039320256561e-13,
                            change_entropy=-0.00024362648639453255,
                            change_entropy_modified=3.240868750253867e-6)
    end

    @trixi_testset "svaerd_kalisch_1d_dingemans_cg" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_dingemans_cg.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.22798490823942433 0.7520004851600044 0.0],
                            linf=[0.03673010870720128 0.12074632168110239 0.0],
                            cons_error=[1.4210854715202004e-13 4.953054817174909e-5 0.0],
                            change_waterheight=-1.4210854715202004e-13,
                            change_entropy=-0.0002425303440531934,
                            change_entropy_modified=2.9729466177741415e-6)
    end

    @trixi_testset "svaerd_kalisch_1d_dingemans_upwind" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_dingemans_upwind.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.2280370308166863 0.7521344942401095 0.0],
                            linf=[0.03673101553812019 0.12116306036094074 0.0],
                            cons_error=[1.1368683772161603e-13 4.871598417571836e-5 0.0],
                            change_waterheight=-1.1368683772161603e-13,
                            change_entropy=-0.00023645232727176335,
                            change_entropy_modified=2.0682244894487667e-6)
    end

    @trixi_testset "svaerd_kalisch_1d_dingemans_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_dingemans_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            N=512,
                            l2=[0.2279601359992612 0.7519295302772953 0.0],
                            linf=[0.03670816991272552 0.12141112705223758 0.0],
                            cons_error=[3.979039320256561e-13 4.937017087502672e-5 0.0],
                            change_waterheight=-3.979039320256561e-13,
                            change_entropy=-0.0002466715038735856,
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
