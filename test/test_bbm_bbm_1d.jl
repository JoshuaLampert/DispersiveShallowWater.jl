module TestBBMBBM1D

using Test
using DispersiveShallowWater

include("test_util.jl")

@testset "BBMBBM1D" begin
  @trixi_testset "bbm_bbm_1d_basic" begin
    trixi_include(@__MODULE__, joinpath(examples_dir(), "bbm_bbm_1d_basic.jl"),
                  tspan = (0.0, 1.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [0.018736231276083863 0.02772945798792677
                    0.016366252270215174 0.0211147869206485
                    2.2432285243130218e-13 1.4779288903810084e-12],
                   atol = 1e-9, rtol = 1e-10)
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [-2.2432285243130218e-13, 1.4779288903810084e-12, 3.729004311026074e-5],
                   atol = 1e-9, rtol = 1e-10)
  end

  @trixi_testset "bbm_bbm_1d_dg" begin
    trixi_include(@__MODULE__, joinpath(examples_dir(), "bbm_bbm_1d_dg.jl"),
                  tspan = (0.0, 1.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [0.14225919400242204 0.04882188307015801
                    0.3421613285081925 0.09882370398079487
                    1.0912164523217173e-15 1.7763568394002505e-15],
                   atol = 1e-9, rtol = 1e-10)
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [1.0912164523217173e-15, 1.7763568394002505e-15, -0.021203943180474027],
                   atol = 1e-8, rtol = 1e-10)
  end

  @trixi_testset "bbm_bbm_1d_relaxation" begin
    trixi_include(@__MODULE__, joinpath(examples_dir(), "bbm_bbm_1d_relaxation.jl"),
                  tspan = (0.0, 1.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [0.018736253012537167 0.02772912154768439
                    0.016366092353111816 0.021114697353509015
                    2.259983767321016e-13 1.4779288903810084e-12],
                   atol = 1e-9, rtol = 1e-10)
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [-2.2461752483080717e-13, 5.684341886080801e-13, 0.0],
                   atol = 1e-9, rtol = 1e-10)
  end
end

end # module
