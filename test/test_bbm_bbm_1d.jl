module TestBBMBBM1D

using Test
using DispersiveShallowWater
using SummationByPartsOperators: integrate

EXAMPLES_DIR = pkgdir(DispersiveShallowWater, "examples")

@testset "BBMBBM1D" begin
  @testset "bbm_bbm_1d_basic" begin
    include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic.jl"))
    analytical_sol = @test_nowarn DispersiveShallowWater.compute_coefficients(initial_condition, solver, tspan[end], equations, mesh)
    l2error_eta =  integrate(u->u^2, sol.u[end][1, :] - analytical_sol[1, :], solver.D) |> sqrt
    @test isapprox(l2error_eta, 0.6761496557998811, atol=1e-7, rtol=1e-8)
    l2error_v =  integrate(u->u^2, sol.u[end][2, :] - analytical_sol[2, :], solver.D) |> sqrt
    @test isapprox(l2error_v, 1.3992447330707167, atol=1e-7, rtol=1e-8)
  end

  @testset "bbm_bbm_1d_dg" begin
    include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_dg.jl"))
    analytical_sol = @test_nowarn DispersiveShallowWater.compute_coefficients(initial_condition, solver, tspan[end], equations, mesh)
    l2error_eta =  integrate(u->u^2, sol.u[end][1, :] - analytical_sol[1, :], solver.D) |> sqrt
    @test isapprox(l2error_eta, 2.029527814597151, atol=1e-7, rtol=1e-8)
    l2error_v =  integrate(u->u^2, sol.u[end][2, :] - analytical_sol[2, :], solver.D) |> sqrt
    @test isapprox(l2error_v, 1.9434053327388578, atol=1e-7, rtol=1e-8)
  end
end

end # module
