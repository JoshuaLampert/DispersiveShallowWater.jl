using Test

@testset "DispersiveShallowWater.jl" begin
  include("test_unit.jl")
  include("test_bbm_bbm_1d.jl")
  include("test_bbm_bbm_variable_bathymetry_1d.jl")
end
