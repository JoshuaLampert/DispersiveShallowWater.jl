module TestBBMBBMVariableBathymetry1D

using Test
using DispersiveShallowWater

include("test_util.jl")

@testset "BBMBBMVariableBathymetry1D" begin
  @trixi_testset "bbm_bbm_1d_variable_bathymetry_basic" begin
    trixi_include(@__MODULE__,
                  joinpath(examples_dir(), "bbm_bbm_variable_bathymetry_1d_basic.jl"),
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

  @trixi_testset "bbm_bbm_1d_variable_bathymetry_relaxation" begin
    trixi_include(@__MODULE__,
                  joinpath(examples_dir(), "bbm_bbm_variable_bathymetry_1d_relaxation.jl"),
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
end

end # module
