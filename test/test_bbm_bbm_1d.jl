module TestBBMBBM1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "bbm_bbm_1d")

@testset "BBMBBM1D" begin
  @trixi_testset "bbm_bbm_1d_basic" begin
    trixi_include(@__MODULE__, joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic.jl"),
                  tspan = (0.0, 1.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [0.018736231276083863 0.02772945798792677
                    0.016366252270215174 0.0211147869206485
                    2.2432285243130218e-13 1.4779288903810084e-12],
                   atol = 500 * eps(), rtol = sqrt(eps()))
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [2.2222469560301384e-13, 1.1368683772161603e-13, 3.729011041286867e-5],
                   atol = 1e-11, rtol = sqrt(eps()))
  end

  @trixi_testset "bbm_bbm_1d_dg" begin
    trixi_include(@__MODULE__, joinpath(EXAMPLES_DIR, "bbm_bbm_1d_dg.jl"),
                  tspan = (0.0, 1.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [0.034635970678256946 0.012225260982110586
                    0.09331575019082416 0.021156308992005712
                    1.0543424823256011e-14 3.552713678800501e-15],
                   atol = 500 * eps(), rtol = sqrt(eps()))
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [1.0543424823256011e-14, -3.552713678800501e-15, -0.043687254604591885],
                   atol = 1e-11, rtol = sqrt(eps()))
  end

  @trixi_testset "bbm_bbm_1d_relaxation" begin
    trixi_include(@__MODULE__, joinpath(EXAMPLES_DIR, "bbm_bbm_1d_relaxation.jl"),
                  tspan = (0.0, 1.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [0.018736253012537167 0.02772912154768439
                    0.016366092353111816 0.021114697353509015
                    2.259983767321016e-13 1.4779288903810084e-12],
                   atol = 500 * eps(), rtol = sqrt(eps()))
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [-2.2461752483080717e-13, 5.684341886080801e-13, 0.0],
                   atol = 1e-11, rtol = sqrt(eps()))
  end
end

end # module
