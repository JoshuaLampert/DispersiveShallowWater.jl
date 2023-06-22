module TestBBMBBMVariableBathymetry1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "bbm_bbm_variable_bathymetry_1d")

@testset "BBMBBMVariableBathymetry1D" begin
  @trixi_testset "bbm_bbm_variable_bathymetry_1d_basic" begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "bbm_bbm_variable_bathymetry_1d_basic.jl"),
                  tspan = (0.0, 1.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [0.04346190082228665 0.11171710004736789 0.0
                    0.04283944204214696 0.0953750884645661 0.0
                    1.0394561372181596e-13 1.1368683772161603e-13 0.0],
                   atol = 500 * eps(), rtol = sqrt(eps()))
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [1.0394561372181596e-13, 1.1368683772161603e-13, 0.06735422396923241],
                   atol = 1e-11, rtol = sqrt(eps()))
  end

  @trixi_testset "bbm_bbm_variable_bathymetry_1d_relaxation" begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "bbm_bbm_variable_bathymetry_1d_relaxation.jl"),
                  tspan = (0.0, 1.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [1.6263014524141195 2.1343189572239996 0.0
                    2.020840307859494 2.6331707343176247 0.0
                    2.0694557179012918e-13 2.789919043027633e-13 0.0],
                   atol = 500 * eps(), rtol = sqrt(eps()))
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [-2.1271873151818e-13, 2.7874259106214573e-13, 0.0],
                   atol = 1e-11, rtol = sqrt(eps()))
  end

  @trixi_testset "bbm_bbm_variable_bathymetry_1d_dg_upwind_relaxation" begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR,
                           "bbm_bbm_variable_bathymetry_1d_dg_upwind_relaxation.jl"),
                  tspan = (0.0, 1.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [1.6263014468802284 2.134318981212618 0.0
                    2.020841127902278 2.6332103195428296 0.0
                    3.2152058793144533e-13 1.0576435560682995e-12 0.0],
                   atol = 500 * eps(), rtol = sqrt(eps()))
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [3.2152058793144533e-13,
                     -1.0576435560682995e-12,
                     -1.4210854715202004e-14],
                   atol = 1e-11, rtol = sqrt(eps()))
  end

  @trixi_testset "bbm_bbm_variable_bathymetry_1d_upwind_relaxation" begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR,
                           "bbm_bbm_variable_bathymetry_1d_upwind_relaxation.jl"),
                  tspan = (0.0, 1.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [1.6263014396164643 2.1343190212669234 0.0
                    2.020820789184524 2.6332272207139242 0.0
                    1.5010215292932116e-13 1.8715064220575783e-13 0.0],
                   atol = 500 * eps(), rtol = sqrt(eps()))
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [-1.5010215292932116e-13, 1.8715064220575783e-13, 0.0],
                   atol = 1e-11, rtol = sqrt(eps()))
  end

  @trixi_testset "bbm_bbm_variable_bathymetry_1d_well_balanced" begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR,
                           "bbm_bbm_variable_bathymetry_1d_well_balanced.jl"),
                  tspan = (0.0, 10.0))
    errs = @view errors(analysis_callback)[:, :, end]
    @test isapprox(errs,
                   [6.313048635750119e-15 3.793418293562131e-14 0.0
                    1.9317880628477724e-14 1.1898289293196742e-13 0.0
                    0.0 5.728935765975126e-15 0.0],
                   atol = 500 * eps(), rtol = sqrt(eps()))
    change_of_invariants = integrals(analysis_callback)[:, end] -
                           integrals(analysis_callback)[:, 1]
    @test isapprox(change_of_invariants,
                   [0.0, -5.728935765975126e-15, 0.0, 6.725397799854067e-15],
                   atol = 1e-11, rtol = sqrt(eps()))
  end
end

end # module