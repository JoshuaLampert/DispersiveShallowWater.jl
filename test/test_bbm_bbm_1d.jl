module TestBBMBBM1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "bbm_bbm_1d")

@testset "BBMBBM1D" begin
    @trixi_testset "bbm_bbm_1d_basic" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.0032744047423353103 0.0072467844128150435],
                            linf=[0.0027330692683324997 0.0037734832582039246],
                            cons_error=[3.61905099297374e-14 5.684341886080802e-14],
                            change_waterheight=3.61905099297374e-14,
                            change_velocity=-5.684341886080802e-14,
                            change_entropy=0.00019552920957721653,
                            atol_ints=1e-10) # in order to make CI pass

        @test_allocations(semi, sol)

        # test upwind operators
        using SummationByPartsOperators: upwind_operators, periodic_derivative_operator
        using SparseArrays: sparse
        using OrdinaryDiffEq: solve
        D1 = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                              accuracy_order = accuracy_order, xmin = mesh.xmin,
                              xmax = mesh.xmax,
                              N = mesh.N)
        D2 = sparse(D1.plus) * sparse(D1.minus)
        solver = Solver(D1, D2)
        semi = Semidiscretization(mesh, equations, initial_condition, solver,
                                  boundary_conditions = boundary_conditions)
        ode = semidiscretize(semi, (0.0, 1.0))
        sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
                    save_everystep = false, callback = callbacks, saveat = saveat)
        atol = 1e-11 # in order to make CI pass
        rtol = 1e-12
        errs = errors(analysis_callback)
        l2 = [0.0024468371786471343 0.004954294231072229]
        l2_measured = errs.l2_error[:, end]
        for (l2_expected, l2_actual) in zip(l2, l2_measured)
            @test isapprox(l2_expected, l2_actual, atol = atol, rtol = rtol)
        end
        linf = [0.002095653257311536 0.0026110374662806635]
        linf_measured = errs.linf_error[:, end]
        for (linf_expected, linf_actual) in zip(linf, linf_measured)
            @test isapprox(linf_expected, linf_actual, atol = atol, rtol = rtol)
        end

        @test_allocations(semi, sol)
    end

    @trixi_testset "bbm_bbm_1d_dg" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_dg.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.034635970678256946 0.012225260982110586],
                            linf=[0.09331575019082416 0.021156308992005712],
                            cons_error=[1.0543424823256011e-14 3.552713678800501e-15],
                            change_waterheight=1.0543424823256011e-14,
                            change_velocity=-3.552713678800501e-15,
                            change_entropy=-0.021843627302246205)

        @test_allocations(semi, sol, allocs=60000)
    end

    @trixi_testset "bbm_bbm_1d_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_relaxation.jl"),
                            tspan=(0.0, 1.0),
                            l2=[0.003273738156062929 0.007245185828684949],
                            linf=[0.0027325632260941646 0.0037724348496581683],
                            cons_error=[2.259983767321016e-13 4.547473508864641e-13],
                            change_waterheight=2.1746572188776938e-13,
                            change_velocity=-4.547473508864641e-13,
                            change_entropy=0.0)

        @test_allocations(semi, sol)
    end

    @trixi_testset "bbm_bbm_1d_manufactured" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_manufactured.jl"),
                            tspan=(0.0, 1.0),
                            l2=[4.365176233405813e-9 6.7151849982388e-10],
                            linf=[6.226559268185383e-9 9.698699621196738e-10],
                            cons_error=[3.873483998828586e-12 2.2986355942039745e-11],
                            change_waterheight=3.873483998828586e-12,
                            change_velocity=2.2986355942039745e-11,
                            change_entropy=17.387441847193436,
                            atol=1e-10,
                            atol_ints=1e-10) # in order to make CI pass

        @test_allocations(semi, sol)
    end

    @trixi_testset "bbm_bbm_1d_basic_reflecting" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_bbm_1d_basic_reflecting.jl"),
                            tspan=(0.0, 1.0),
                            l2=[1.4463385687002652e-6 2.5981071881535818e-8],
                            linf=[2.6974309994542978e-6 4.159023703209641e-8],
                            cons_error=[4.272316442782926e-11 0.5469460931577768],
                            change_waterheight=4.272316442782926e-11,
                            change_velocity=0.5469460931577768,
                            change_entropy=130.69415963528348)

        @test_allocations(semi, sol)

        # test upwind operators
        using SummationByPartsOperators: upwind_operators, Mattsson2017
        using SparseArrays: sparse
        using OrdinaryDiffEq: solve
        D1 = upwind_operators(Mattsson2017; derivative_order = 1,
                              accuracy_order = accuracy_order, xmin = mesh.xmin,
                              xmax = mesh.xmax,
                              N = mesh.N)
        D2 = sparse(D1.plus) * sparse(D1.minus)
        solver = Solver(D1, D2)
        semi = Semidiscretization(mesh, equations, initial_condition, solver,
                                  boundary_conditions = boundary_conditions,
                                  source_terms = source_terms)
        ode = semidiscretize(semi, (0.0, 1.0))
        sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
                    save_everystep = false, callback = callbacks, saveat = saveat)
        atol = 1e-10 # in order to make CI pass
        rtol = 1e-12
        errs = errors(analysis_callback)
        l2 = [6.46559852118817e-6 2.2682195103606648e-8]
        l2_measured = errs.l2_error[:, end]
        for (l2_expected, l2_actual) in zip(l2, l2_measured)
            @test isapprox(l2_expected, l2_actual, atol = atol, rtol = rtol)
        end
        linf = [0.00015506981854507274 8.639888084832625e-8]
        linf_measured = errs.linf_error[:, end]
        for (linf_expected, linf_actual) in zip(linf, linf_measured)
            @test isapprox(linf_expected, linf_actual, atol = atol, rtol = rtol)
        end

        @test_allocations(semi, sol)
    end
end

end # module
