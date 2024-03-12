module TestAqua

using Aqua
using ExplicitImports: check_no_implicit_imports, check_no_stale_explicit_imports
using Test
using DispersiveShallowWater

@testset "Aqua.jl" begin
    Aqua.test_all(DispersiveShallowWater,
                  ambiguities = false)
    @test isnothing(check_no_implicit_imports(DispersiveShallowWater))
    @test isnothing(check_no_stale_explicit_imports(DispersiveShallowWater))
end

end #module 
