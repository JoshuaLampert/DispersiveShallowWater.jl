module TestAqua

using Aqua
using Test
using DispersiveShallowWater

@testset "Aqua.jl" begin
    Aqua.test_all(DispersiveShallowWater,
                  ambiguities = false)
end

end #module 
