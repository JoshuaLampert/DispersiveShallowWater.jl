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
                   [0.01889143455446012 0.02942294364283494
                    0.017059188018341255 0.021405665120430406
                    2.2426934041972823e-13 1.2505552149377763e-12],
                   atol = 1e-9, rtol = 1e-10)
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [-2.2426934041972823e-13, 1.2505552149377763e-12, 0.014085559731029207],
                   atol = 1e-9, rtol = 1e-10)
  end

  @trixi_testset "bbm_bbm_1d_dg" begin
    trixi_include(@__MODULE__, joinpath(examples_dir(), "bbm_bbm_1d_dg.jl"),
                  tspan = (0.0, 1.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [0.14214505895648474 0.04880262293404891
                    0.3419093520999896 0.09948660703480527
                    2.7819374137216107e-15 3.552713678800501e-15],
                   atol = 1e-9, rtol = 1e-10)
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [2.7819374137216107e-15, -3.552713678800501e-15, 0.0011275899092524355],
                   atol = 1e-8, rtol = 1e-10)
  end

  @trixi_testset "bbm_bbm_1d_relaxation" begin
    trixi_include(@__MODULE__, joinpath(examples_dir(), "bbm_bbm_1d_relaxation.jl"),
                  tspan = (0.0, 1.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [0.018885586625155475 0.0292345575077295
                    0.016976265438477856 0.021374041350339823
                    2.246215224604027e-13 1.4779288903810084e-12],
                   atol = 1e-9, rtol = 1e-10)
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [-2.246215224604027e-13, 1.4779288903810084e-12, 0.0],
                   atol = 1e-9, rtol = 1e-10)
  end
end

end # module
