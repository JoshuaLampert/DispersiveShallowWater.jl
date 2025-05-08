@testitem "Aqua.jl" setup=[Setup] begin
    using Aqua
    using ExplicitImports: check_no_implicit_imports, check_no_stale_explicit_imports
    using JET

    # Aqua.jl
    Aqua.test_all(DispersiveShallowWater,
                  ambiguities = false)

    # ExplicitImports.jl
    @test isnothing(check_no_implicit_imports(DispersiveShallowWater))
    @test isnothing(check_no_stale_explicit_imports(DispersiveShallowWater))

    # JET.jl
    # With the default settings as of 2025-05-08, JET.jl
    # reports issues originating from RecipesBase.jl.
    # The only way to ignore them seems to be something like the
    # following hack inspired by the discussion in
    # https://github.com/aviatesk/JET.jl/issues/570
    struct IgnoreRecipesBase end
    function JET.match_module(::IgnoreRecipesBase,
                              @nospecialize(report::JET.InferenceErrorReport))
        s = "MethodInstance for RecipesBase.apply_recipe"
        any(report.vst) do vst
            occursin(s, string(vst))
        end
    end
    test_package(DispersiveShallowWater;
                 target_defined_modules = true,
                 ignored_modules = (IgnoreRecipesBase(),))
end
