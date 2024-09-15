module TestBBMBBM1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "bbm_1d")

@testset "BBMEquation1D" begin
    @trixi_testset "bbm_1d_basic" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_basic.jl"),
                            tspan=(0.0, 100.0),
                            l2=[0.0002514795321988484],
                            linf=[9.756568022900591e-5],
                            cons_error=[7.105427357601002e-15],
                            change_waterheight=-7.105427357601002e-15,
                            change_entropy_modified=-4.4274395039067826e-7,
                            change_hamiltonian=-5.33011109249415e-7)

        @test_allocations(semi, sol, allocs=5_000)

        # test upwind operators
        using SummationByPartsOperators: upwind_operators, periodic_derivative_operator
        using SparseArrays: sparse
        using OrdinaryDiffEq: solve
        D1 = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                              accuracy_order = accuracy_order, xmin = mesh.xmin,
                              xmax = mesh.xmax,
                              N = mesh.N)
        D2 = sparse(D1.minus) * sparse(D1.plus)
        solver = Solver(D1, D2)
        semi = Semidiscretization(mesh, equations, initial_condition, solver,
                                  boundary_conditions = boundary_conditions)
        ode = semidiscretize(semi, (0.0, 100.0))
        sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
                    save_everystep = false, callback = callbacks, saveat = saveat)
        atol = 1e-12
        rtol = 1e-12
        errs = errors(analysis_callback)
        l2 = [0.00018397011159873]
        l2_measured = errs.l2_error[:, end]
        for (l2_expected, l2_actual) in zip(l2, l2_measured)
            @test isapprox(l2_expected, l2_actual, atol = atol, rtol = rtol)
        end
        linf = [7.190886191621448e-5]
        linf_measured = errs.linf_error[:, end]
        for (linf_expected, linf_actual) in zip(linf, linf_measured)
            @test isapprox(linf_expected, linf_actual, atol = atol, rtol = rtol)
        end

        @test_allocations(semi, sol, allocs=5_000)
    end

    @trixi_testset "bbm_1d_basic with split_form = false" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_basic.jl"),
                            tspan=(0.0, 100.0),
                            split_form=false,
                            l2=[0.00023679720962102973],
                            linf=[9.291484780826753e-5],
                            cons_error=[0.0],
                            change_waterheight=0.0,
                            change_entropy_modified=-4.4114923292148944e-7,
                            change_hamiltonian=-5.313053028643822e-7)

        @test_allocations(semi, sol, allocs=5_000)

        # test upwind operators
        using SummationByPartsOperators: upwind_operators, periodic_derivative_operator
        using SparseArrays: sparse
        using OrdinaryDiffEq: solve
        D1 = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                              accuracy_order = accuracy_order, xmin = mesh.xmin,
                              xmax = mesh.xmax,
                              N = mesh.N)
        D2 = sparse(D1.minus) * sparse(D1.plus)
        solver = Solver(D1, D2)
        semi = Semidiscretization(mesh, equations, initial_condition, solver,
                                  boundary_conditions = boundary_conditions)
        ode = semidiscretize(semi, (0.0, 100.0))
        sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
                    save_everystep = false, callback = callbacks, saveat = saveat)
        atol = 1e-12
        rtol = 1e-12
        errs = errors(analysis_callback)
        l2 = [0.00016247570402673777]
        l2_measured = errs.l2_error[:, end]
        for (l2_expected, l2_actual) in zip(l2, l2_measured)
            @test isapprox(l2_expected, l2_actual, atol = atol, rtol = rtol)
        end
        linf = [6.495267031642049e-5]
        linf_measured = errs.linf_error[:, end]
        for (linf_expected, linf_actual) in zip(linf, linf_measured)
            @test isapprox(linf_expected, linf_actual, atol = atol, rtol = rtol)
        end

        @test_allocations(semi, sol, allocs=5_000)
    end

    @trixi_testset "bbm_1d_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_relaxation.jl"),
                            tspan=(0.0, 100.0),
                            l2=[0.0002510482044505331],
                            linf=[9.745861147153478e-5],
                            cons_error=[5.329070518200751e-15],
                            change_waterheight=5.329070518200751e-15,
                            change_entropy_modified=0.0,
                            change_hamiltonian=-1.7181174261082788e-9)

        @test_allocations(semi, sol, allocs=5_000)
    end

    @trixi_testset "bbm_1d_hamiltonian_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_hamiltonian_relaxation.jl"),
                            tspan=(0.0, 100.0),
                            l2=[0.00023637369561147064],
                            linf=[9.280778692660752e-5],
                            cons_error=[1.7763568394002505e-14],
                            change_waterheight=-1.7763568394002505e-14,
                            change_entropy_modified=1.6049757078917537e-9,
                            change_hamiltonian=0.0)

        @test_allocations(semi, sol, allocs=5_000)
    end

    @trixi_testset "bbm_1d_fourier" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_fourier.jl"),
                            tspan=(0.0, 100.0),
                            l2=[6.01265210291245e-7],
                            linf=[2.3226949846799627e-7],
                            cons_error=[1.1546319456101628e-14],
                            change_waterheight=-1.1546319456101628e-14,
                            change_entropy_modified=-4.4266038123907947e-7,
                            change_hamiltonian=-5.311924822226644e-7)

        @test_allocations(semi, sol, allocs=5_000)
    end

    @trixi_testset "bbm_1d_manufactured" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_manufactured.jl"),
                            tspan=(0.0, 1.0),
                            l2=[7.887878295265557e-9],
                            linf=[9.041914550422803e-9],
                            cons_error=[5.148179246187295e-13],
                            change_waterheight=5.148179246187295e-13)

        @test_allocations(semi, sol, allocs=5_000)

        # test upwind operators
        using SummationByPartsOperators: upwind_operators, periodic_derivative_operator
        using SparseArrays: sparse
        using OrdinaryDiffEq: solve
        D1 = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                              accuracy_order = accuracy_order, xmin = mesh.xmin,
                              xmax = mesh.xmax,
                              N = mesh.N)
        D2 = sparse(D1.minus) * sparse(D1.plus)
        solver = Solver(D1, D2)
        semi = Semidiscretization(mesh, equations, initial_condition, solver,
                                  boundary_conditions = boundary_conditions,
                                  source_terms = source_terms)
        ode = semidiscretize(semi, (0.0, 1.0))
        sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
                    save_everystep = false, callback = callbacks, saveat = saveat)
        atol = 1e-12
        rtol = 1e-12
        errs = errors(analysis_callback)
        l2 = [9.784591180473132e-8]
        l2_measured = errs.l2_error[:, end]
        for (l2_expected, l2_actual) in zip(l2, l2_measured)
            @test isapprox(l2_expected, l2_actual, atol = atol, rtol = rtol)
        end
        linf = [9.915897769552373e-8]
        linf_measured = errs.linf_error[:, end]
        for (linf_expected, linf_actual) in zip(linf, linf_measured)
            @test isapprox(linf_expected, linf_actual, atol = atol, rtol = rtol)
        end

        @test_allocations(semi, sol, allocs=5_000)
    end
end

end # module
