module TestSerreGreenNaghdi1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "serre_green_naghdi_1d")

@testset "SerreGreenNaghdiEquations1D" begin
    @trixi_testset "serre_green_naghdi_soliton.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton.jl"),
                            tspan=(0.0, 0.1),
                            l2=[9.994998669268741e-7, 1.4703973445698635e-6, 0.0],
                            linf=[6.5496216650196e-7, 1.027617322124641e-6, 0.0],
                            cons_error=[0.0, 8.174581012099225e-10, 0.0],
                            change_waterheight=0.0,
                            change_entropy_modified=-3.1093350116861984e-11)

        @test_allocations(semi, sol, 550_000)
    end

    @trixi_testset "serre_green_naghdi_soliton.jl with bathymetry_mild_slope" begin
        # same values as serre_green_naghdi_soliton.jl but with more allocations
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton.jl"),
                            equations = SerreGreenNaghdiEquations1D(bathymetry_mild_slope; gravity_constant = 9.81),
                            tspan=(0.0, 0.1),
                            l2=[9.994998669268741e-7, 1.4703973445698635e-6, 0.0],
                            linf=[6.5496216650196e-7, 1.027617322124641e-6, 0.0],
                            cons_error=[0.0, 8.174581012099225e-10, 0.0],
                            change_waterheight=0.0,
                            change_entropy_modified=-3.1093350116861984e-11)

        @test_allocations(semi, sol, 800_000)
    end

    @trixi_testset "serre_green_naghdi_soliton.jl with bathymetry_variable" begin
        # same values as serre_green_naghdi_soliton.jl but with more allocations
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton.jl"),
                            equations = SerreGreenNaghdiEquations1D(bathymetry_variable; gravity_constant = 9.81),
                            tspan=(0.0, 0.1),
                            l2=[9.994998669268741e-7, 1.4703973445698635e-6, 0.0],
                            linf=[6.5496216650196e-7, 1.027617322124641e-6, 0.0],
                            cons_error=[0.0, 8.174581012099225e-10, 0.0],
                            change_waterheight=0.0,
                            change_entropy_modified=-3.1093350116861984e-11)

        @test_allocations(semi, sol, 800_000)
    end

    @trixi_testset "serre_green_naghdi_soliton_fourier.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton_fourier.jl"),
                            tspan=(0.0, 0.1),
                            l2=[8.252225014546995e-8, 6.724994492548714e-7, 0.0],
                            linf=[2.672093302180656e-8, 9.642725156897014e-8, 0.0],
                            cons_error=[2.842170943040401e-14, 4.627409566637652e-13, 0.0],
                            change_waterheight=2.842170943040401e-14,
                            change_entropy_modified=-3.097966327914037e-11)

        @test_allocations(semi, sol, allocs=450_000)
    end

    @trixi_testset "serre_green_naghdi_soliton_fourier.jl with bathymetry_mild_slope" begin
        # same values as serre_green_naghdi_soliton_fourier.jl but with more allocations
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton_fourier.jl"),
                            equations = SerreGreenNaghdiEquations1D(bathymetry_mild_slope; gravity_constant = 9.81),
                            tspan=(0.0, 0.1),
                            l2=[8.252225014546995e-8, 6.724994492548714e-7, 0.0],
                            linf=[2.672093302180656e-8, 9.642725156897014e-8, 0.0],
                            cons_error=[2.842170943040401e-14, 4.627409566637652e-13, 0.0],
                            change_waterheight=2.842170943040401e-14,
                            change_entropy_modified=-3.097966327914037e-11)

        @test_allocations(semi, sol, allocs=850_000)
    end

    @trixi_testset "serre_green_naghdi_soliton_fourier.jl with bathymetry_variable" begin
        # same values as serre_green_naghdi_soliton_fourier.jl but with more allocations
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton_fourier.jl"),
                            equations = SerreGreenNaghdiEquations1D(bathymetry_variable; gravity_constant = 9.81),
                            tspan=(0.0, 0.1),
                            l2=[8.252225014546995e-8, 6.724994492548714e-7, 0.0],
                            linf=[2.672093302180656e-8, 9.642725156897014e-8, 0.0],
                            cons_error=[2.842170943040401e-14, 4.627409566637652e-13, 0.0],
                            change_waterheight=2.842170943040401e-14,
                            change_entropy_modified=-3.097966327914037e-11)

        @test_allocations(semi, sol, allocs=850_000)
    end

    @trixi_testset "serre_green_naghdi_soliton_upwind.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton_upwind.jl"),
                            tspan=(0.0, 0.1),
                            l2=[1.4876412924488654e-6, 5.9888097605442856e-6, 0.0],
                            linf=[1.0863034516361836e-6, 4.105927902009476e-6, 0.0],
                            cons_error=[4.263256414560601e-14, 4.483030568991353e-8, 0.0],
                            change_waterheight=4.263256414560601e-14,
                            change_entropy_modified=-3.1036506698001176e-11)

        @test_allocations(semi, sol, allocs=500_000)
    end

    @trixi_testset "serre_green_naghdi_soliton_upwind.jl with bathymetry_mild_slope" begin
        # same as serre_green_naghdi_soliton_upwind.jl but with more allocations
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton_upwind.jl"),
                            tspan=(0.0, 0.1),
                            equations = SerreGreenNaghdiEquations1D(bathymetry_mild_slope; gravity_constant = 9.81),
                            l2=[1.4876412924488654e-6, 5.9888097605442856e-6, 0.0],
                            linf=[1.0863034516361836e-6, 4.105927902009476e-6, 0.0],
                            cons_error=[4.263256414560601e-14, 4.483030568991353e-8, 0.0],
                            change_waterheight=4.263256414560601e-14,
                            change_entropy_modified=-3.1036506698001176e-11)

        @test_allocations(semi, sol, allocs=750_000)
    end

    @trixi_testset "serre_green_naghdi_soliton_upwind.jl with bathymetry_variable" begin
        # same as serre_green_naghdi_soliton_upwind.jl but with more allocations
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton_upwind.jl"),
                            tspan=(0.0, 0.1),
                            equations = SerreGreenNaghdiEquations1D(bathymetry_variable; gravity_constant = 9.81),
                            l2=[1.4876412924488654e-6, 5.9888097605442856e-6, 0.0],
                            linf=[1.0863034516361836e-6, 4.105927902009476e-6, 0.0],
                            cons_error=[4.263256414560601e-14, 4.483030568991353e-8, 0.0],
                            change_waterheight=4.263256414560601e-14,
                            change_entropy_modified=-3.1036506698001176e-11)

        @test_allocations(semi, sol, allocs=750_000)
    end

    @trixi_testset "serre_green_naghdi_soliton_relaxation.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton_relaxation.jl"),
                            tspan=(0.0, 0.1),
                            l2=[8.252225169608892e-8, 6.724994488577288e-7, 0.0],
                            linf=[2.6716495016287922e-8, 9.642466235713909e-8, 0.0],
                            cons_error=[2.842170943040401e-14, 4.649614027130156e-13, 0.0],
                            change_waterheight=2.842170943040401e-14,
                            change_entropy_modified=0.0)

        @test_allocations(semi, sol, allocs=450_000)
    end

    @trixi_testset "serre_green_naghdi_soliton_relaxation.jl with bathymetry_mild_slope" begin
        # same as serre_green_naghdi_soliton_relaxation.jl but with more allocations
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton_relaxation.jl"),
                            tspan=(0.0, 0.1),
                            equations = SerreGreenNaghdiEquations1D(bathymetry_mild_slope; gravity_constant = 9.81),
                            l2=[8.252225169608892e-8, 6.724994488577288e-7, 0.0],
                            linf=[2.6716495016287922e-8, 9.642466235713909e-8, 0.0],
                            cons_error=[2.842170943040401e-14, 4.649614027130156e-13, 0.0],
                            change_waterheight=2.842170943040401e-14,
                            change_entropy_modified=0.0)

        @test_allocations(semi, sol, allocs=850_000)
    end

    @trixi_testset "serre_green_naghdi_soliton_relaxation.jl with bathymetry_variable" begin
        # same as serre_green_naghdi_soliton_relaxation.jl but with more allocations
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_soliton_relaxation.jl"),
                            tspan=(0.0, 0.1),
                            equations = SerreGreenNaghdiEquations1D(bathymetry_variable; gravity_constant = 9.81),
                            l2=[8.252225169608892e-8, 6.724994488577288e-7, 0.0],
                            linf=[2.6716495016287922e-8, 9.642466235713909e-8, 0.0],
                            cons_error=[2.842170943040401e-14, 4.649614027130156e-13, 0.0],
                            change_waterheight=2.842170943040401e-14,
                            change_entropy_modified=0.0)

        @test_allocations(semi, sol, allocs=850_000)
    end

    @trixi_testset "serre_green_naghdi_dingemans.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_dingemans.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.2463295492331802, 0.805357513755897, 0.0],
                            linf=[0.036265772921736605, 0.11763218152350845, 0.0],
                            cons_error=[5.684341886080802e-14, 3.509473432346808e-5, 0.0],
                            change_waterheight=5.684341886080802e-14,
                            change_entropy=1.9489684518703143e-5,
                            change_entropy_modified=-1.2177679309388623e-8)

        @test_allocations(semi, sol, allocs=750_000)
    end

    @trixi_testset "serre_green_naghdi_dingemans.jl with bathymetry_mild_slope" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "serre_green_naghdi_dingemans.jl"),
                            tspan=(0.0, 1.0),
                            equations = SerreGreenNaghdiEquations1D(bathymetry_mild_slope; gravity_constant = 9.81),
                            l2=[0.2463295492331802, 0.805357513755897, 0.0],
                            linf=[0.036265772921736605, 0.11763218152350845, 0.0],
                            cons_error=[5.684341886080802e-14, 3.509473432346624e-5, 0.0],
                            change_waterheight=5.684341886080802e-14,
                            change_entropy=1.9489684518703143e-5,
                            change_entropy_modified=-1.2177679309388623e-8)

        @test_allocations(semi, sol, allocs=750_000)
    end
end

end # module
