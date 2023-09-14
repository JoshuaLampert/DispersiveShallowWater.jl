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
                            l2=[0.0024561065856105136 0.005467472679563469 0.0],
                            linf=[0.0019322426802252934 0.003270790905389731 0.0],
                            cons_error=[1.0394561372181596e-13 3.410605131648481e-13 0.0],
                            change_waterheight=1.0394561372181596e-13,
                            change_velocity=-3.410605131648481e-13,
                            change_entropy=0.0006138716858004045)
    end

    @trixi_testset "bbm_bbm_variable_bathymetry_1d_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[1.234497154595864 2.0325352963894225 0.0],
                            linf=[1.856609837506367 2.537322089418664 0.0],
                            cons_error=[2.0694557179012918e-13 2.789919043027633e-13 0.0],
                            change_waterheight=-2.1271873151818e-13,
                            change_velocity=2.7874259106214573e-13,
                            change_entropy=0.0)
    end

    @trixi_testset "bbm_bbm_variable_bathymetry_1d_dg_upwind_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_dg_upwind_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[1.2344971565355631 2.0325352892395023 0.0],
                            linf=[1.8563719900575752 2.537410881290326 0.0],
                            cons_error=[3.1796787425264483e-13 1.080125572316959e-14 0.0],
                            change_waterheight=3.1796787425264483e-13,
                            change_velocity=-1.0584519372081047e-12,
                            change_entropy=1.4210854715202004e-14)
    end

    @trixi_testset "bbm_bbm_variable_bathymetry_1d_upwind_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_upwind_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[1.2344971597115844 2.032535283544955 0.0],
                            linf=[1.8566997938473593 2.5374108392574604 0.0],
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

    @trixi_testset "bbm_bbm_variable_bathymetry_1d_dingemans" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_dingemans.jl"),
                            tspan=(0.0, 1.0),
                            N=512,
                            l2=[0.23120795192990629 0.7544344847900466 0.0],
                            linf=[0.03707736413915552 0.12038327721384329 0.0],
                            cons_error=[1.4210854715202004e-13 3.1478183774857893e-15 0.0],
                            change_waterheight=-1.4210854715202004e-13,
                            change_velocity=-3.1478183774857893e-15,
                            change_entropy=-1.1468728189356625e-9)
    end
end

end # module
