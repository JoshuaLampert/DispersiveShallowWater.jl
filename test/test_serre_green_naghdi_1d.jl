module TestSerreGreenNaghdi1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "serre_green_naghdi_1d")

@testset "SerreGreenNaghdiEquations1D" begin
    @trixi_testset "serre_green_naghdi_soliton.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton.jl"),
                            tspan=(0.0, 0.1),
                            l2=[9.994998669268741e-7, 1.4703955789342244e-6],
                            linf=[6.5496216650196e-7, 1.027615853854691e-6],
                            cons_error=[0.0, 8.174581012099225e-10],
                            change_waterheight=0.0,
                            change_entropy_modified=-3.1093350116861984e-11,
                            atol=1e-9) # in order to make CI pass

        @test_allocations(semi, sol, allocs=550_000)
    end

    @trixi_testset "serre_green_naghdi_soliton_fourier.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton_fourier.jl"),
                            tspan=(0.0, 0.1),
                            l2=[8.252225014546995e-8, 6.724994492548714e-7],
                            linf=[2.6719845003242426e-8, 9.642725156897014e-8],
                            cons_error=[2.842170943040401e-14, 4.627409566637652e-13],
                            change_waterheight=2.842170943040401e-14,
                            change_entropy_modified=-3.097966327914037e-11,
                            atol=1e-9) # in order to make CI pass

        @test_allocations(semi, sol, allocs=550_000)
    end

    @trixi_testset "serre_green_naghdi_soliton_upwind.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton_upwind.jl"),
                            tspan=(0.0, 0.1),
                            l2=[1.4876412924488654e-6, 5.988810995995645e-6],
                            linf=[1.0863034516361836e-6, 4.105929242048667e-6],
                            cons_error=[4.263256414560601e-14, 4.483030568991353e-8],
                            change_waterheight=4.263256414560601e-14,
                            change_entropy_modified=-3.1036506698001176e-11,
                            atol=1e-9) # in order to make CI pass

        @test_allocations(semi, sol, allocs=550_000)
    end

    @trixi_testset "serre_green_naghdi_soliton_relaxation.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton_relaxation.jl"),
                            tspan=(0.0, 0.1),
                            l2=[8.252225169608892e-8, 6.724994488577288e-7],
                            linf=[2.6722121404532118e-8, 9.64274304918189e-8],
                            cons_error=[2.842170943040401e-14, 4.649614027130156e-13],
                            change_waterheight=2.842170943040401e-14,
                            change_entropy_modified=0.0,
                            atol=1e-9) # in order to make CI pass

        @test_allocations(semi, sol, allocs=550_000)
    end
end

end # module
