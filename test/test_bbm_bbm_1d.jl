module TestBBMBBM1D

using Test
using DispersiveShallowWater

include("test_util.jl")

@testset "BBMBBM1D" begin
  @trixi_testset "bbm_bbm_1d_basic" begin
    using SummationByPartsOperators: integrate

    trixi_include(@__MODULE__, joinpath(examples_dir(), "bbm_bbm_1d_basic.jl"),
                  tspan = (0.0, 1.0))
    analytical_sol = @test_nowarn DispersiveShallowWater.compute_coefficients(initial_condition,
                                                                              solver,
                                                                              tspan[end],
                                                                              equations,
                                                                              mesh)
    l2error_eta = integrate(u -> u^2, sol.u[end][1, :] - analytical_sol[1, :], solver.D) |>
                  sqrt
    @test isapprox(l2error_eta, 0.01989008182686233, atol = 1e-9, rtol = 1e-10)
    l2error_v = integrate(u -> u^2, sol.u[end][2, :] - analytical_sol[2, :], solver.D) |>
                sqrt
    @test isapprox(l2error_v, 0.035239757311485494, atol = 1e-9, rtol = 1e-10)
  end

  @trixi_testset "bbm_bbm_1d_dg" begin
    using SummationByPartsOperators: integrate

    trixi_include(@__MODULE__, joinpath(examples_dir(), "bbm_bbm_1d_dg.jl"),
                  tspan = (0.0, 1.0))
    analytical_sol = @test_nowarn DispersiveShallowWater.compute_coefficients(initial_condition,
                                                                              solver,
                                                                              tspan[end],
                                                                              equations,
                                                                              mesh)
    l2error_eta = integrate(u -> u^2, sol.u[end][1, :] - analytical_sol[1, :], solver.D) |>
                  sqrt
    @test isapprox(l2error_eta, 0.14214505895648474, atol = 1e-9, rtol = 1e-10)
    l2error_v = integrate(u -> u^2, sol.u[end][2, :] - analytical_sol[2, :], solver.D) |>
                sqrt
    @test isapprox(l2error_v, 0.04880262293404891, atol = 1e-9, rtol = 1e-10)
  end
end

end # module
