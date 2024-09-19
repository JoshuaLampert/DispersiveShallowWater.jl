@testitem "Aqua.jl" setup=[Setup] begin
    using Aqua
    using ExplicitImports: check_no_implicit_imports, check_no_stale_explicit_imports
    Aqua.test_all(DispersiveShallowWater,
                  ambiguities = false)
    @test isnothing(check_no_implicit_imports(DispersiveShallowWater))
    @test isnothing(check_no_stale_explicit_imports(DispersiveShallowWater))
end
