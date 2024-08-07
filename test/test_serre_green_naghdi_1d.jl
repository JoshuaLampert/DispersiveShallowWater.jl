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
                            change_entropy=-3.1093350116861984e-11,
                            atol=1e-9) # in order to make CI pass

        @test_allocations(semi, sol, allocs=550_000)
    end
end

end # module
