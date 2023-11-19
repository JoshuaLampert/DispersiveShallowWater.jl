module TestBBMBBM1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "bbm_bbm_1d")

@testset "BBMBBM1D" begin
    @trixi_testset "bbm_bbm_1d_basic" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.0032744047432744098 0.007246784414780245],
                            linf=[0.0027330692690079594 0.0037734832590992085],
                            cons_error=[2.234687408354354e-13 5.684341886080801e-13],
                            change_waterheight=2.2222469560301384e-13,
                            change_velocity=-5.684341886080801e-13,
                            change_entropy=0.00019552914864107152,
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
                            change_entropy=-0.021843627302246205)
    end

    @trixi_testset "bbm_bbm_1d_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.003273738156062929 0.007245185828684949],
                            linf=[0.0027325632260941646 0.0037724348496581683],
                            cons_error=[2.259983767321016e-13 4.547473508864641e-13],
                            change_waterheight=2.1746572188776938e-13,
                            change_velocity=-4.547473508864641e-13,
                            change_entropy=0.0)
    end

    @trixi_testset "bbm_bbm_1d_manufactured" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_manufactured.jl"),
                            tspan=(0.0, 1.0),
                            l2=[4.365176233405813e-9 6.7151849982388e-10],
                            linf=[6.226559268185383e-9 9.698699621196738e-10],
                            cons_error=[3.873483998828586e-12 2.2986355942039745e-11],
                            change_waterheight=3.873483998828586e-12,
                            change_velocity=2.2986355942039745e-11,
                            change_entropy=17.387441847193436,
                            atol=1e-10,
                            atol_ints=1e-10) # in order to make CI pass
    end
end

end # module
