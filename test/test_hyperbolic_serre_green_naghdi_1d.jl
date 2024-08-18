module TestSerreGreenNaghdi1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "hyperbolic_serre_green_naghdi_1d")

@testset "HyperbolicSerreGreenNaghdiEquations1D" begin
    @trixi_testset "hyperbolic_serre_green_naghdi_soliton.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "hyperbolic_serre_green_naghdi_soliton.jl"),
                            tspan=(0.0, 0.1),
                            l2=[0.0007038714283663042, 0.006508261273448058, 0.0, 0.024517136865274798, 0.002141907410685252],
                            linf=[0.0005088046605401519, 0.0036954890877776703, 0.0, 0.015022422297545818, 0.0013290414555349184],
                            cons_error=[2.7000623958883807e-13, 0.00013389587974454997, 0.0, 0.005963937086921899, 4.502801745331908e-5],
                            change_entropy_modified=-2.3374946067633573e-7)

        @test_allocations(semi, sol, 1_000)
    end

    @trixi_testset "hyperbolic_serre_green_naghdi_soliton_relaxation.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "hyperbolic_serre_green_naghdi_soliton_relaxation.jl"),
                            tspan=(0.0, 0.1),
                            l2=[0.0007041797674417557, 0.006510737539763134, 0.0, 0.024517447804525746, 0.002141928791106223],
                            linf=[0.0005090662088376163, 0.0036987746989370907, 0.0, 0.01502004552088677, 0.0013289272946777064],
                            cons_error=[3.126388037344441e-13, 0.00013409338344283483, 0.0, 0.005963706457799891, 4.504121848469822e-5],
                            change_entropy_modified=-5.684341886080802e-14)

        @test_allocations(semi, sol, 1_000)
    end
end

end # module
