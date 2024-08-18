using Test

@testset "DispersiveShallowWater.jl" begin
    include("test_aqua.jl")
    include("test_unit.jl")
    include("test_visualization.jl")
    include("test_bbm_bbm_1d.jl")
    include("test_bbm_bbm_variable_bathymetry_1d.jl")
    include("test_svaerd_kalisch_1d.jl")
    include("test_serre_green_naghdi_1d.jl")
    include("test_hyperbolic_serre_green_naghdi_1d.jl")
end
