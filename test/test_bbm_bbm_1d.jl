@testsnippet BBMBBMEquation1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "bbm_bbm_1d")
end

@testitem "bbm_bbm_1d_basic" setup=[Setup, BBMBBMEquation1D, AdditionalImports] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.0032743987817220976 0.007246750130235789 0.0],
                        linf=[0.0027330553477815656 0.0037734958705577526 0.0],
                        cons_error=[2.7030772060575454e-14 1.4210854715202004e-13 0.0],
                        change_waterheight=-2.7030772060575454e-14,
                        change_velocity=-1.4210854715202004e-13,
                        change_entropy=0.00023829378642403753,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=10_000)

    # test upwind operators
    D1 = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                          accuracy_order = accuracy_order, xmin = mesh.xmin,
                          xmax = mesh.xmax,
                          N = mesh.N)
    D2 = sparse(D1.plus) * sparse(D1.minus)
    solver = Solver(D1, D2)
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic.jl"),
                        tspan=(0.0, 1.0),
                        solver=solver,
                        l2=[0.0026452454364373377 0.005181138085234744 0.0],
                        linf=[0.002240297682717385 0.003138940449190386 0.0],
                        cons_error=[6.490606933236491e-14 3.694822225952521e-13 0.0],
                        change_waterheight=6.490606933236491e-14,
                        change_velocity=3.694822225952521e-13,
                        change_entropy=0.0002383181188179151,
                        atol_ints=1e-11) # to make CI pass

    @test_allocations(semi, sol, allocs=10_000)

    # test PeriodicRationalDerivativeOperator
    D1 = periodic_derivative_operator(1, accuracy_order, xmin(mesh), xmax(mesh),
                                      nnodes(mesh))
    D2 = D1^2
    solver = Solver(D1, D2)
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic.jl"),
                        tspan=(0.0, 1.0),
                        solver=solver,
                        l2=[0.0016951648276611925 0.0034269307395771303 0.0],
                        linf=[0.001453478820216958 0.001805665634286413 0.0],
                        cons_error=[2.055466838899258e-14 0.0 0.0],
                        change_waterheight=2.055466838899258e-14,
                        change_velocity=0.0,
                        change_entropy=0.00023833941781958856,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_basic with bathymetry_variable" setup=[Setup, BBMBBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic.jl"),
                        bathymetry_type=bathymetry_variable,
                        tspan=(0.0, 1.0),
                        l2=[0.002475254552148715 0.005510109246665293 0.0],
                        linf=[0.0019397546359500861 0.0032957358113563373 0.0],
                        cons_error=[1.0394561372181596e-13 3.410605131648481e-13 0.0],
                        change_waterheight=1.0394561372181596e-13,
                        change_velocity=-3.410605131648481e-13,
                        change_entropy=0.0006197172569955001,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_dg" setup=[Setup, BBMBBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_dg.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.005820370895289719 0.0012483148286926665 0.0],
                        linf=[0.012036829933978588 0.0022436858082413025 0.0],
                        cons_error=[2.2555754226055478e-14 7.815970093361102e-14 0.0],
                        change_waterheight=2.2555754226055478e-14,
                        change_velocity=-7.815970093361102e-14,
                        change_entropy=-0.001165974183408025)

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_dg with bathymetry_variable" setup=[Setup, BBMBBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_dg.jl"),
                        bathymetry_type=bathymetry_variable,
                        tspan=(0.0, 1.0),
                        l2=[0.009197531778215928 0.0015521194273725176 0.0],
                        linf=[0.012034113656467671 0.0033862062614615773 0.0],
                        cons_error=[3.364296388121304e-14 9.947598300641403e-14 0.0],
                        change_waterheight=3.364296388121304e-14,
                        change_velocity=-9.947598300641403e-14,
                        change_entropy=-0.000791230828923517)

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_fourier" setup=[Setup, BBMBBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_fourier.jl"),
                        tspan=(0.0, 1.0),
                        l2=[7.291521968885259e-7 7.702204594455217e-7 0.0],
                        linf=[4.3579646391567195e-7 4.9442969896063e-7 0.0],
                        cons_error=[3.9206314028412105e-14 4.547473508864641e-13 0.0],
                        change_waterheight=3.9206314028412105e-14,
                        change_velocity=-4.547473508864641e-13,
                        change_entropy=0.0002383147650562023,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_fourier with bathymetry_variable" setup=[Setup, BBMBBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_fourier.jl"),
                        bathymetry_type=bathymetry_variable,
                        tspan=(0.0, 1.0),
                        l2=[7.291481123937497e-7 7.70214005540768e-7 0.0],
                        linf=[4.357912374297612e-7 4.94414834406598e-7 0.0],
                        cons_error=[6.232593470137382e-13 4.547473508864641e-13 0.0],
                        change_waterheight=6.232593470137382e-13,
                        change_velocity=-4.547473508864641e-13,
                        change_entropy=0.0002383154298968293,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_relaxation" setup=[Setup, BBMBBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_relaxation.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.0032736306642711533 0.007244906041050075 0.0],
                        linf=[0.002732472407469544 0.0037722741328174436 0.0],
                        cons_error=[2.47001411123548e-14 2.842170943040401e-14 0.0],
                        change_waterheight=-2.47001411123548e-14,
                        change_velocity=2.842170943040401e-14,
                        change_entropy=-2.9558577807620168e-12)

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_relaxation with bathymetry_variable" setup=[Setup, BBMBBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_relaxation.jl"),
                        bathymetry_type=bathymetry_variable,
                        tspan=(0.0, 1.0),
                        l2=[0.0024773815236814826 0.005515225594228142 0.0],
                        linf=[0.001941460860594768 0.0032992043482167333 0.0],
                        cons_error=[2.864588178487867e-14 5.684341886080802e-14 0.0],
                        change_waterheight=2.864588178487867e-14,
                        change_velocity=-5.684341886080802e-14,
                        change_entropy=-2.2737367544323206e-13)

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_upwind_relaxation" setup=[Setup, BBMBBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_upwind_relaxation.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.002644536846506582 0.005179311517589997 0.0],
                        linf=[0.002239772357125247 0.0031376618560692293 0.0],
                        cons_error=[3.1799390576478594e-14 3.979039320256561e-13 0.0],
                        change_waterheight=3.1799390576478594e-14,
                        change_velocity=3.979039320256561e-13,
                        change_entropy=6.821210263296962e-13)

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_upwind_relaxation with bathymetry_variable" setup=[
    Setup,
    BBMBBMEquation1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_upwind_relaxation.jl"),
                        bathymetry_type=bathymetry_variable,
                        tspan=(0.0, 1.0),
                        l2=[0.00264453684734065 0.0051793115198710705 0.0],
                        linf=[0.0022397723578102546 0.0031376618574547877 0.0],
                        cons_error=[4.234505741492232e-15 4.547473508864641e-13 0.0],
                        change_waterheight=4.234505741492232e-15,
                        change_velocity=4.547473508864641e-13,
                        change_entropy=-2.2737367544323206e-12)

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_manufactured" setup=[Setup, BBMBBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_manufactured.jl"),
                        tspan=(0.0, 1.0),
                        l2=[6.866821527172215e-9 3.4469275186089985e-9 0.0],
                        linf=[1.1689290868588387e-8 5.928427526669111e-9 0.0],
                        cons_error=[2.1126156379835948e-12 8.243777181889383e-12 0.0],
                        change_waterheight=-2.1126156379835948e-12,
                        change_velocity=-8.243777181889383e-12,
                        change_entropy=17.817012269328853,
                        atol=1e-10,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_manufactured with bathymetry_flat" setup=[Setup, BBMBBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_manufactured.jl"),
                        bathymetry_type=bathymetry_flat,
                        tspan=(0.0, 1.0),
                        l2=[5.564654145503227e-9 8.058808980392836e-10 0.0],
                        linf=[7.905728516277577e-9 1.1420695500419242e-9 0.0],
                        cons_error=[5.767395241940143e-12 4.184121331471304e-12 0.0],
                        change_waterheight=-5.767395241940143e-12,
                        change_velocity=4.184121331471304e-12,
                        change_entropy=17.38744182699684,
                        atol=1e-11,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_basic_reflecting" setup=[Setup, BBMBBMEquation1D, AdditionalImports] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic_reflecting.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.00022223640611354037 7.626931465504548e-9 0.0],
                        linf=[0.0002947437933196184 1.7383174233387422e-8 0.0],
                        cons_error=[9.701410286113879e-10 0.5469460962472683 0.0],
                        change_waterheight=9.701410286113879e-10,
                        change_velocity=0.5469460962472683,
                        change_entropy=132.10935771957952,
                        atol=1e-10,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=1_500)

    # test upwind operators
    D1 = upwind_operators(Mattsson2017; derivative_order = 1,
                          accuracy_order = accuracy_order, xmin = mesh.xmin,
                          xmax = mesh.xmax,
                          N = mesh.N)
    D2 = sparse(D1.plus) * sparse(D1.minus)
    solver = Solver(D1, D2)
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic_reflecting.jl"),
                        tspan=(0.0, 1.0),
                        solver=solver,
                        l2=[0.002345799818043513 3.254313503127441e-8 0.0],
                        linf=[0.05625950533062252 6.815531732318192e-7 0.0],
                        cons_error=[1.6607871307809518e-9 0.5469460993745239 0.0],
                        change_waterheight=-1.6607871307809518e-9,
                        change_velocity=0.5469460993745239,
                        change_entropy=132.10938489083918,
                        atol=1e-7,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=2_000)
end

@testitem "bbm_bbm_1d_basic_reflecting with bathymetry_flat" setup=[
    Setup,
    BBMBBMEquation1D,
    AdditionalImports
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic_reflecting.jl"),
                        bathymetry_type=bathymetry_flat,
                        tspan=(0.0, 1.0),
                        l2=[7.585650541014755e-6 5.3346100010644425e-9 0.0],
                        linf=[1.016759621386143e-5 9.171166714949663e-9 0.0],
                        cons_error=[1.13200113340443e-10 0.5469460971120386 0.0],
                        change_waterheight=1.13200113340443e-10,
                        change_velocity=0.5469460971120386,
                        change_entropy=132.04866881559724,
                        atol=1e-10,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=1_500)

    # test upwind operators
    D1 = upwind_operators(Mattsson2017; derivative_order = 1,
                          accuracy_order = accuracy_order, xmin = mesh.xmin,
                          xmax = mesh.xmax,
                          N = mesh.N)
    D2 = sparse(D1.plus) * sparse(D1.minus)
    solver = Solver(D1, D2)
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic_reflecting.jl"),
                        bathymetry_type=bathymetry_flat,
                        tspan=(0.0, 1.0),
                        solver=solver,
                        l2=[9.791289242253764e-5 6.608214058228884e-9 0.0],
                        linf=[0.0023482494155127043 8.684315878915161e-8 0.0],
                        cons_error=[1.4778607518007766e-10 0.5469460970051805 0.0],
                        change_waterheight=-1.4778607518007766e-10,
                        change_velocity=0.5469460970051805,
                        change_entropy=132.04866880504562,
                        atol=1e-8,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=2_000)
end

@testitem "bbm_bbm_1d_dingemans" setup=[Setup, BBMBBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_dingemans.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.22292157345226027 0.7504924411607958 0.0],
                        linf=[0.03584030574058168 0.1202292994661 0.0],
                        cons_error=[2.6129272356900657e-17 2.1141942363467336e-16 0.0],
                        change_waterheight=-2.6129272356900657e-17,
                        change_velocity=2.1141942363467336e-16,
                        change_entropy=3.4175334942543323e-7)

    @test_allocations(semi, sol, allocs=10_000)
end

@testitem "bbm_bbm_1d_well_balanced" setup=[Setup, BBMBBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_well_balanced.jl"),
                        tspan=(0.0, 1.0),
                        l2=[5.625847040853183e-15 2.7315317493074205e-15 0.0],
                        linf=[5.995204332975845e-15 1.12866211315047e-14 0.0],
                        cons_error=[8.881784197001252e-16 1.6076981933671723e-15 0.0],
                        change_waterheight=-8.881784197001252e-16,
                        change_velocity=-1.6076981933671723e-15,
                        change_entropy=-3.1086244689504383e-15)

    @test_allocations(semi, sol, allocs=10_000)
end
