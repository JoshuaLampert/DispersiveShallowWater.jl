@testsnippet SerreGreenNaghdiEquations1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "serre_green_naghdi_1d")
end

@testitem "serre_green_naghdi_soliton.jl" setup=[Setup, SerreGreenNaghdiEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton.jl"),
                        tspan=(0.0, 0.1),
                        l2=[9.994998669268741e-7, 1.4703973445698635e-6, 0.0],
                        linf=[6.5496216650196e-7, 1.027617322124641e-6, 0.0],
                        cons_error=[0.0, 8.174581012099225e-10, 0.0],
                        change_waterheight=0.0,
                        change_entropy_modified=-3.1093350116861984e-11)

    @test_allocations(semi, sol, allocs=550_000)
end

@testitem "serre_green_naghdi_soliton.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same values as serre_green_naghdi_soliton.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        tspan=(0.0, 0.1),
                        l2=[9.994998669268741e-7, 1.4703973445698635e-6, 0.0],
                        linf=[6.5496216650196e-7, 1.027617322124641e-6, 0.0],
                        cons_error=[0.0, 8.174581012099225e-10, 0.0],
                        change_waterheight=0.0,
                        change_entropy_modified=-3.1093350116861984e-11)

    @test_allocations(semi, sol, allocs=800_000)
end

@testitem "serre_green_naghdi_soliton.jl with bathymetry_variable" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same values as serre_green_naghdi_soliton.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton.jl"),
                        bathymetry_type=bathymetry_variable,
                        tspan=(0.0, 0.1),
                        l2=[9.994998669268741e-7, 1.4703973445698635e-6, 0.0],
                        linf=[6.5496216650196e-7, 1.027617322124641e-6, 0.0],
                        cons_error=[0.0, 8.174581012099225e-10, 0.0],
                        change_waterheight=0.0,
                        change_entropy_modified=-3.1093350116861984e-11)

    @test_allocations(semi, sol, allocs=800_000)
end

@testitem "serre_green_naghdi_soliton_fourier.jl" setup=[Setup, SerreGreenNaghdiEquations1D] begin
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

@testitem "serre_green_naghdi_soliton_fourier.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same values as serre_green_naghdi_soliton_fourier.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_fourier.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        tspan=(0.0, 0.1),
                        l2=[8.252225014546995e-8, 6.724994492548714e-7, 0.0],
                        linf=[2.672093302180656e-8, 9.642725156897014e-8, 0.0],
                        cons_error=[2.842170943040401e-14, 4.627409566637652e-13, 0.0],
                        change_waterheight=2.842170943040401e-14,
                        change_entropy_modified=-3.097966327914037e-11)

    @test_allocations(semi, sol, allocs=850_000)
end

@testitem "serre_green_naghdi_soliton_fourier.jl with bathymetry_variable" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same values as serre_green_naghdi_soliton_fourier.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_fourier.jl"),
                        bathymetry_type=bathymetry_variable,
                        tspan=(0.0, 0.1),
                        l2=[8.252225014546995e-8, 6.724994492548714e-7, 0.0],
                        linf=[2.672093302180656e-8, 9.642725156897014e-8, 0.0],
                        cons_error=[2.842170943040401e-14, 4.627409566637652e-13, 0.0],
                        change_waterheight=2.842170943040401e-14,
                        change_entropy_modified=-3.097966327914037e-11)

    @test_allocations(semi, sol, allocs=850_000)
end

@testitem "serre_green_naghdi_soliton_upwind.jl" setup=[Setup, SerreGreenNaghdiEquations1D] begin
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

@testitem "serre_green_naghdi_soliton_upwind.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same as serre_green_naghdi_soliton_upwind.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_upwind.jl"),
                        tspan=(0.0, 0.1),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[1.4876412924488654e-6, 5.9888097605442856e-6, 0.0],
                        linf=[1.0863034516361836e-6, 4.105927902009476e-6, 0.0],
                        cons_error=[4.263256414560601e-14, 4.483030568991353e-8, 0.0],
                        change_waterheight=4.263256414560601e-14,
                        change_entropy_modified=-3.1036506698001176e-11)

    @test_allocations(semi, sol, allocs=750_000)
end

@testitem "serre_green_naghdi_soliton_upwind.jl with bathymetry_variable" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same as serre_green_naghdi_soliton_upwind.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_upwind.jl"),
                        tspan=(0.0, 0.1),
                        bathymetry_type=bathymetry_variable,
                        l2=[1.4876412924488654e-6, 5.9888097605442856e-6, 0.0],
                        linf=[1.0863034516361836e-6, 4.105927902009476e-6, 0.0],
                        cons_error=[4.263256414560601e-14, 4.483030568991353e-8, 0.0],
                        change_waterheight=4.263256414560601e-14,
                        change_entropy_modified=-3.1036506698001176e-11)

    @test_allocations(semi, sol, allocs=750_000)
end

@testitem "serre_green_naghdi_soliton_relaxation.jl" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
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

@testitem "serre_green_naghdi_soliton_relaxation.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same as serre_green_naghdi_soliton_relaxation.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_relaxation.jl"),
                        tspan=(0.0, 0.1),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[8.252225169608892e-8, 6.724994488577288e-7, 0.0],
                        linf=[2.6716495016287922e-8, 9.642466235713909e-8, 0.0],
                        cons_error=[2.842170943040401e-14, 4.649614027130156e-13, 0.0],
                        change_waterheight=2.842170943040401e-14,
                        change_entropy_modified=0.0)

    @test_allocations(semi, sol, allocs=850_000)
end

@testitem "serre_green_naghdi_soliton_relaxation.jl with bathymetry_variable" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same as serre_green_naghdi_soliton_relaxation.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_relaxation.jl"),
                        tspan=(0.0, 0.1),
                        bathymetry_type=bathymetry_variable,
                        l2=[8.252225169608892e-8, 6.724994488577288e-7, 0.0],
                        linf=[2.6716495016287922e-8, 9.642466235713909e-8, 0.0],
                        cons_error=[2.842170943040401e-14, 4.649614027130156e-13, 0.0],
                        change_waterheight=2.842170943040401e-14,
                        change_entropy_modified=0.0)

    @test_allocations(semi, sol, allocs=850_000)
end

@testitem "serre_green_naghdi_well_balanced.jl" setup=[Setup, SerreGreenNaghdiEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_well_balanced.jl"),
                        tspan=(0.0, 2.0),
                        l2=[0, 0, 0],
                        linf=[0, 0, 0],
                        cons_error=[0, 0, 0],
                        change_waterheight=0.0,
                        change_momentum=0.0,
                        change_entropy_modified=0.0,
                        lake_at_rest=0.0)

    @test_allocations(semi, sol, allocs=750_000)
end

@testitem "serre_green_naghdi_well_balanced.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_well_balanced.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        tspan=(0.0, 2.0),
                        l2=[0, 0, 0],
                        linf=[0, 0, 0],
                        cons_error=[0, 0, 0],
                        change_waterheight=0.0,
                        change_momentum=0.0,
                        change_entropy_modified=0.0,
                        lake_at_rest=0.0)

    @test_allocations(semi, sol, allocs=750_000)
end

@testitem "serre_green_naghdi_dingemans.jl" setup=[Setup, SerreGreenNaghdiEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_dingemans.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.22632930215585131, 0.7400070292134782, 0.0],
                        linf=[0.036351214376643126, 0.11899056101300992, 0.0],
                        cons_error=[1.4210854715202004e-13, 3.194346928167053e-5, 0.0],
                        change_waterheight=-1.4210854715202004e-13,
                        change_entropy=2.282635693973134e-5,
                        change_entropy_modified=-9.135646905633621e-9)

    @test_allocations(semi, sol, allocs=750_000)
end

@testitem "serre_green_naghdi_dingemans.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_dingemans.jl"),
                        tspan=(0.0, 1.0),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[0.22632930215585131, 0.7400070292134782, 0.0],
                        linf=[0.036351214376643126, 0.11899056101300992, 0.0],
                        cons_error=[1.4210854715202004e-13, 3.194346928167053e-5, 0.0],
                        change_waterheight=-1.4210854715202004e-13,
                        change_entropy=2.282635693973134e-5,
                        change_entropy_modified=-9.135646905633621e-9)

    @test_allocations(semi, sol, allocs=750_000)
end

@testitem "serre_green_naghdi_conservation.jl" setup=[Setup, SerreGreenNaghdiEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_conservation.jl"),
                        l2=[1.3655498085989206, 2.3967486930606716, 0.0],
                        linf=[1.001076318001934, 0.8052527556023067, 0.0],
                        cons_error=[0.0, 0.0002674927404067162, 0.0],
                        change_entropy=-0.05841897226287074,
                        change_entropy_modified=0.059273551933074486,
                        atol_ints=2e-8, # to make CI pass
                        atol=2e-8) # to make CI pass

    @test_allocations(semi, sol, allocs=900_000)
end

@testitem "serre_green_naghdi_conservation.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_conservation.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[1.3655493671985637, 2.3967828251339003, 0.0],
                        linf=[1.001075913983051, 0.8052680970114169, 0.0],
                        cons_error=[1.1368683772161603e-13, 0.00026407261543415217, 0.0],
                        change_entropy=-0.058352273553509804,
                        change_entropy_modified=0.05927340849780194,
                        atol_ints=2e-8) # to make CI pass

    @test_allocations(semi, sol, allocs=900_000)
end
