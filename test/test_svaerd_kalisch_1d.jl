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
                            l2=[0.22648616606810756 0.7325813437573053 0.0],
                            linf=[0.03645879748792702 0.11849678361651249 0.0],
                            cons_error=[3.979039320256561e-13 7.670877062280329e-5 0.0],
                            change_waterheight=-3.979039320256561e-13,
                            change_momentum=-2.109259079896564e-9,
                            change_entropy=-0.0005910026394531087,
                            change_entropy_modified=3.20516357987799e-6)
    end

    @trixi_testset "svaerd_kalisch_1d_dingemans_upwind" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_dingemans_upwind.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.22656501380746077 0.7327645752975466 0.0],
                            linf=[0.03648151661791921 0.1182504804050214 0.0],
                            cons_error=[1.1368683772161603e-13 7.59953335162671e-5 0.0],
                            change_waterheight=-1.1368683772161603e-13,
                            change_momentum=-2.28489619585881e-9,
                            change_entropy=-0.0005840425510541536,
                            change_entropy_modified=2.108449052684591e-6)
    end

    @trixi_testset "svaerd_kalisch_1d_dingemans_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_dingemans_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            N=512,
                            l2=[0.22648560566941997 0.7325790781930053 0.0],
                            linf=[0.03645867682442583 0.1184963632507331 0.0],
                            cons_error=[3.979039320256561e-13 7.670526860626775e-5 0.0],
                            change_waterheight=-3.979039320256561e-13,
                            change_entropy=-0.0005940058100577517,
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
