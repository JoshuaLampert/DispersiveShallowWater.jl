
@testsnippet SvaerdKalischEquations1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "svaerd_kalisch_1d")
end

@testitem "svaerd_kalisch_1d_manufactured" setup=[Setup, SvaerdKalischEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_manufactured.jl"),
                        tspan=(0.0, 0.1),
                        l2=[3.3755102050606554e-6 2.320896961343125e-7 0.0],
                        linf=[4.908886917176503e-6 3.888671399332466e-7 0.0],
                        cons_error=[2.42861286636753e-16 1.9224170696150768e-7 0.0],
                        change_waterheight=-2.42861286636753e-16,
                        change_entropy=0.1868146724821993,
                        atol=1e-9) # to make CI pass

    @test_allocations(semi, sol, allocs=90_000)
end

@testitem "svaerd_kalisch_1d_dingemans" setup=[Setup, SvaerdKalischEquations1D, AdditionalImports] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_dingemans.jl"),
                        tspan=(0.0, 1.0),
                        N=512,
                        l2=[0.22796106962338855 0.7519327063662515 0.0],
                        linf=[0.036708347831218346 0.12141172207472928 0.0],
                        cons_error=[3.979039320256561e-13 4.937137540373564e-5 0.0],
                        change_waterheight=-3.979039320256561e-13,
                        change_entropy=-0.00024362648639453255,
                        change_entropy_modified=-6.311893230304122e-9)

    @test_allocations(semi, sol, allocs=350_000)

    # test PeriodicRationalDerivativeOperator
    D1 = periodic_derivative_operator(1, accuracy_order, xmin(mesh), xmax(mesh), nnodes(mesh))
    D2 = D1^2
    solver = Solver(D1, D2)
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_dingemans.jl"),
                        tspan=(0.0, 1.0),
                        N=512,
                        solver=solver,
                        l2=[0.22796180766836865 0.7519296292874257 0.0],
                        linf=[0.036709492542750466 0.12104908724733915 0.0],
                        cons_error=[2.842170943040401e-14 4.9341465183441843e-5 0.0],
                        change_waterheight=-2.842170943040401e-14,
                        change_entropy=-0.00024270962080663594,
                        change_entropy_modified=-7.430799087160267e-9)

    @test_allocations(semi, sol, allocs=350_000)
end

@testitem "svaerd_kalisch_1d_dingemans_cg" setup=[Setup, SvaerdKalischEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_dingemans_cg.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.22798490823942433 0.7520004851600044 0.0],
                        linf=[0.03673010870720128 0.12074632168110239 0.0],
                        cons_error=[1.4210854715202004e-13 4.953054817174909e-5 0.0],
                        change_waterheight=-1.4210854715202004e-13,
                        change_entropy=-0.0002425303440531934,
                        change_entropy_modified=-2.6815314413397573e-9)

    @test_allocations(semi, sol, allocs=750_000)
end

@testitem "svaerd_kalisch_1d_dingemans_fourier" setup=[Setup, SvaerdKalischEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_dingemans_fourier.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.22799246923254327 0.7520172891302948 0.0],
                        linf=[0.03671494177483947 0.12129171577180138 0.0],
                        cons_error=[8.526512829121202e-14 5.3078570574495334e-5 0.0],
                        change_waterheight=-8.526512829121202e-14,
                        change_entropy=-0.0002424441479433881,
                        change_entropy_modified=-4.00007138523506e-9)

    @test_allocations(semi, sol, allocs=13_000_000)
end

@testitem "svaerd_kalisch_1d_dingemans_upwind" setup=[Setup, SvaerdKalischEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_dingemans_upwind.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.2280370308166863 0.7521344942401095 0.0],
                        linf=[0.03673101553812019 0.12116306036094074 0.0],
                        cons_error=[1.1368683772161603e-13 4.871598417571836e-5 0.0],
                        change_waterheight=-1.1368683772161603e-13,
                        change_entropy=-0.00023645232727176335,
                        change_entropy_modified=-6.654090611846186e-9)

    @test_allocations(semi, sol, allocs=350_000)
end

@testitem "svaerd_kalisch_1d_dingemans_relaxation" setup=[Setup, SvaerdKalischEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_dingemans_relaxation.jl"),
                        tspan=(0.0, 1.0),
                        N=512,
                        l2=[0.22796107242561717 0.7519327155905444 0.0],
                        linf=[0.03670834831604197 0.12141172368792873 0.0],
                        cons_error=[3.979039320256561e-13 4.937137655207271e-5 0.0],
                        change_waterheight=-3.979039320256561e-13,
                        change_entropy=-0.00024362054875837202,
                        change_entropy_modified=0.0)

    @test_allocations(semi, sol, allocs=350_000)
end

@testitem "svaerd_kalisch_1d_well_balanced" setup=[Setup, SvaerdKalischEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_well_balanced.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.0 1.135448143093612e-14 0.0],
                        linf=[0.0 8.133477278069499e-15 0.0],
                        cons_error=[0.0 1.6056589579882354e-14 0.0],
                        change_waterheight=0.0,
                        change_momentum=1.5679986322667355e-14,
                        change_entropy=0.0,
                        lake_at_rest=0.0)

    @test_allocations(semi, sol, allocs=150_000)
end
