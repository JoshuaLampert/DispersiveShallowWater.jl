module TestUnit

using Test
using DispersiveShallowWater
using Plots

@testset "Visualization" begin
    trixi_include(@__MODULE__, default_example(), tspan = (0.0, 1.0))
    @test_nowarn plot(semi => sol)
    @test_nowarn plot!(semi => sol, plot_initial = true)
    @test_nowarn plot(semi, sol, conversion = prim2cons, plot_bathymetry = false)
    @test_nowarn plot(semi => sol, 0.0)
    @test_nowarn plot(semi, sol, 0.0, conversion = prim2cons)
    @test_nowarn plot(analysis_callback)
    @test_nowarn plot(analysis_callback, what = (:errors,))
    @test_nowarn plot(analysis_callback, what = (:integrals, :errors))
end

end # module
