using Test

@testset "DispersiveShallowWater.jl" begin
  include("test_unit.jl")
  include("test_bbm_bbm_1d.jl")
end
