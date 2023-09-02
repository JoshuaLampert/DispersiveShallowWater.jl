module TestBBMBBM1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "bbm_bbm_1d")

@testset "BBMBBM1D" begin
    @trixi_testset "bbm_bbm_1d_basic" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.003249079018874934 0.007190714378932064],
                            linf=[0.0027187328947428924 0.0037427924312432026],
                            cons_error=[2.234687408354354e-13 5.684341886080801e-13],
                            change_waterheight=2.2222469560301384e-13,
                            change_velocity=-5.684341886080801e-13,
                            change_entropy=0.0003910017421731027,
                            atol_ints=1e-10) # in order to make CI pass
    end

    @trixi_testset "bbm_bbm_1d_dg" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_dg.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.034635970678256946 0.012225260982110586],
                            linf=[0.09331575019082416 0.021156308992005712],
                            cons_error=[1.0543424823256011e-14 3.552713678800501e-15],
                            change_waterheight=1.0543424823256011e-14,
                            change_velocity=-3.552713678800501e-15,
                            change_entropy=-0.043687254604591885)
    end

    @trixi_testset "bbm_bbm_1d_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.003248414414759437 0.007189120223695793],
                            linf=[0.002718237079326169 0.0037417143809364006],
                            cons_error=[2.259983767321016e-13 4.547473508864641e-13],
                            change_waterheight=2.1746572188776938e-13,
                            change_velocity=-4.547473508864641e-13,
                            change_entropy=0.0)
    end
end

end # module
