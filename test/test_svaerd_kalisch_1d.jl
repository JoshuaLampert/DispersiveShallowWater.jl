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
                            l2=[0.2264997669219955 0.7326153404586632 0.0],
                            linf=[0.03641032763299723 0.11738950764868172 0.0],
                            cons_error=[3.979039320256561e-13 7.648176553013466e-5 0.0],
                            change_waterheight=-3.979039320256561e-13,
                            change_momentum=-2.0323272786892943e-9,
                            change_entropy=-0.0005894453857990811,
                            change_entropy_modified=3.0048763619561214e-6)
    end

    @trixi_testset "svaerd_kalisch_1d_dingemans_upwind" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_dingemans_upwind.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.2265797299117127 0.7327963925248041 0.0],
                            linf=[0.036426119662156475 0.1178108264598402 0.0],
                            cons_error=[1.1368683772161603e-13 7.540163948846978e-5 0.0],
                            change_waterheight=-1.1368683772161603e-13,
                            change_momentum=-2.1784861692353275e-9,
                            change_entropy=-0.0005782276325589919,
                            change_entropy_modified=1.7219122128153685e-6)
    end

    @trixi_testset "svaerd_kalisch_1d_dingemans_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "svaerd_kalisch_1d_dingemans_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            N=512,
                            l2=[0.22649926749071012 0.7326132878628044 0.0],
                            linf=[0.03641021393801758 0.11738893742931189 0.0],
                            cons_error=[3.979039320256561e-13 7.647827842231048e-5 0.0],
                            change_waterheight=-3.979039320256561e-13,
                            change_entropy=-0.0005922598786582967,
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
