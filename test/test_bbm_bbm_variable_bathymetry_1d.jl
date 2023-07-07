module TestBBMBBMVariableBathymetry1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "bbm_bbm_variable_bathymetry_1d")

@testset "BBMBBMVariableBathymetry1D" begin
    @trixi_testset "bbm_bbm_variable_bathymetry_1d_basic" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_basic.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.04346190082228665 0.11171710004736789 0.0],
                            linf=[0.04283944204214696 0.0953750884645661 0.0],
                            cons_error=[1.0394561372181596e-13 3.410605131648481e-13 0.0],
                            change_waterheight=1.0394561372181596e-13,
                            change_velocity=-3.410605131648481e-13,
                            change_entropy=0.06735422396923241)
    end

    @trixi_testset "bbm_bbm_variable_bathymetry_1d_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[1.6263014524141195 2.1343189572239996 0.0],
                            linf=[2.020840307859494 2.6331707343176247 0.0],
                            cons_error=[2.0694557179012918e-13 2.789919043027633e-13 0.0],
                            change_waterheight=-2.1271873151818e-13,
                            change_velocity=2.7874259106214573e-13,
                            change_entropy=0.0)
    end

    @trixi_testset "bbm_bbm_variable_bathymetry_1d_dg_upwind_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_dg_upwind_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[1.6263014468792105 2.134318981217363 0.0],
                            linf=[2.020841121130176 2.6332103195755705 0.0],
                            cons_error=[3.1796787425264483e-13 1.0584519372081047e-12 0.0],
                            change_waterheight=3.1796787425264483e-13,
                            change_velocity=-1.0584519372081047e-12,
                            change_entropy=1.4210854715202004e-14)
    end

    @trixi_testset "bbm_bbm_variable_bathymetry_1d_upwind_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_upwind_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[1.626301439652868 2.134319022870026 0.0],
                            linf=[2.0208207914664107 2.633227224286964 0.0],
                            cons_error=[7.194245199571014e-14 1.2388735870505485e-12 0.0],
                            change_waterheight=-1.5765166949677223e-13,
                            change_velocity=1.8791999206735355e-13,
                            change_entropy=2.1316282072803006e-14)
    end

    @trixi_testset "bbm_bbm_variable_bathymetry_1d_well_balanced" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_well_balanced.jl"),
                            tspan=(0.0, 1.0),
                            l2=[6.313048635750119e-15 3.793418293562131e-14 0.0],
                            linf=[1.9317880628477724e-14 1.1898289293196742e-13 0.0],
                            cons_error=[0.0 5.728935765975126e-15 0.0],
                            change_waterheight=0.0,
                            change_velocity=-5.728935765975126e-15,
                            change_entropy=0.0,
                            lake_at_rest=6.725397799854067e-15)
    end
end

end # module
