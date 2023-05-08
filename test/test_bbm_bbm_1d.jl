module TestBBMBBM1D

using Test
using DispersiveShallowWater
using SummationByPartsOperators: integrate

EXAMPLES_DIR = joinpath(pathof(DispersiveShallowWater) |> dirname |> dirname, "examples")

# Start with a clean environment: remove output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive=true)

@testset "BBMBBM1D" begin
  @testset "bbm_bbm_1d_basic" begin
    include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic.jl"))
    analytical_sol = @test_nowarn DispersiveShallowWater.compute_coefficients(initial_condition, solver, tspan[end], equations, mesh)
    l2error_eta =  integrate(u->u^2, sol.u[end][1, :] - analytical_sol[1, :], solver.D) |> sqrt
    @test isapprox(l2error_eta, 0.5155722604645978, atol=500*eps(), rtol=sqrt(eps()))
    l2error_v =  integrate(u->u^2, sol.u[end][2, :] - analytical_sol[2, :], solver.D) |> sqrt
    @test isapprox(l2error_v, 0.48900191864693715, atol=500*eps(), rtol=sqrt(eps()))
  end

  @testset "bbm_bbm_1d_dg" begin
    include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_dg.jl"))
    analytical_sol = @test_nowarn DispersiveShallowWater.compute_coefficients(initial_condition, solver, tspan[end], equations, mesh)
    l2error_eta =  integrate(u->u^2, sol.u[end][1, :] - analytical_sol[1, :], solver.D) |> sqrt
    @test isapprox(l2error_eta, 2.029527814597151, atol=500*eps(), rtol=sqrt(eps()))
    l2error_v =  integrate(u->u^2, sol.u[end][2, :] - analytical_sol[2, :], solver.D) |> sqrt
    @test isapprox(l2error_v, 1.9434053327388578, atol=500*eps(), rtol=sqrt(eps()))
  end
end

# Clean up afterwards: delete output directory
@test_nowarn rm(outdir, recursive=true)

end # module
