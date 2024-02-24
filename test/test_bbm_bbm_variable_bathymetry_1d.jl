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
                            l2=[0.002475254552148715 0.005510109246665293 0.0],
                            linf=[0.0019397546359500861 0.0032957358113563373 0.0],
                            cons_error=[1.0394561372181596e-13 3.410605131648481e-13 0.0],
                            change_waterheight=1.0394561372181596e-13,
                            change_velocity=-3.410605131648481e-13,
                            change_entropy=0.0006197172569955001,
                            atol_ints=1e-10) # in order to make CI pass
    end

    @trixi_testset "bbm_bbm_variable_bathymetry_1d_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.001695855504730616 0.0034287779832081213 0.0],
                            linf=[0.0014540429957374812 0.001806888012062302 0.0],
                            cons_error=[2.0694557179012918e-13 2.789919043027633e-13 0.0],
                            change_waterheight=-2.1271873151818e-13,
                            change_velocity=2.7874259106214573e-13,
                            change_entropy=0.0)
    end

    @trixi_testset "bbm_bbm_variable_bathymetry_1d_dg_upwind_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_dg_upwind_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.004960390675877118 0.00830158298558452 0.0],
                            linf=[0.011473564092181476 0.016650286103717254 0.0],
                            cons_error=[3.1796787425264483e-13 1.080125572316959e-14 0.0],
                            change_waterheight=3.1796787425264483e-13,
                            change_velocity=-1.0584519372081047e-12,
                            change_entropy=1.4210854715202004e-14)
    end

    @trixi_testset "bbm_bbm_variable_bathymetry_1d_upwind_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_upwind_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.002644536848587577 0.005179311522968244 0.0],
                            linf=[0.0022397723587049834 0.0031376618593412786 0.0],
                            cons_error=[7.194245199571014e-14 1.2388735870505485e-12 0.0],
                            change_waterheight=-1.5765166949677223e-13,
                            change_velocity=1.8791999206735355e-13,
                            change_entropy=2.1316282072803006e-14,
                            atol=1e-11) # in order to make CI pass)
    end

    @trixi_testset "bbm_bbm_variable_bathymetry_1d_manufactured" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "bbm_bbm_variable_bathymetry_1d_manufactured.jl"),
                            tspan=(0.0, 1.0),
                            l2=[6.867287425380051e-9 3.446178245128195e-9 0.0],
                            linf=[1.1667709465257303e-8 5.917459633408839e-9 0.0],
                            cons_error=[1.359075785939412e-11 3.8711139735371144e-13 0.0],
                            change_waterheight=-1.359075785939412e-11,
                            change_velocity=-3.8711139735371144e-13,
                            change_entropy=17.81701226932122,
                            atol=1e-10,
                            atol_ints=1e-11) # in order to make CI pass
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
                            l2=[0.2229215734522602 0.750492441160796 0.0],
                            linf=[0.03584030574058169 0.1202292994661 0.0],
                            cons_error=[1.4210854715202004e-13 3.1478183774857893e-15 0.0],
                            change_waterheight=-1.4210854715202004e-13,
                            change_velocity=-3.1478183774857893e-15,
                            change_entropy=3.417533493699221e-7)
    end
end

end # module
