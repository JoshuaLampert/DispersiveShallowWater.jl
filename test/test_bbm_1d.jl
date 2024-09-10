module TestBBMBBM1D

using Test
using DispersiveShallowWater

include("test_util.jl")

EXAMPLES_DIR = joinpath(examples_dir(), "bbm_1d")

@testset "BBM1D" begin
    @trixi_testset "bbm_1d_basic" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_basic.jl"),
                            tspan=(0.0, 100.0),
                            l2=[0.0002514795321988484],
                            linf=[9.756568022900591e-5],
                            cons_error=[7.105427357601002e-15],
                            change_waterheight=-7.105427357601002e-15,
                            change_entropy_modified=6.355536918301041e-6,
                            change_invariant_cubic=-3.198066679033218e-6)

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

    @trixi_testset "bbm_1d_relaxation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_relaxation.jl"),
                            tspan=(0.0, 100.0),
                            l2=[0.00026259713456399205],
                            linf=[0.00010095858421543813],
                            cons_error=[6.217248937900877e-15],
                            change_waterheight=6.217248937900877e-15,
                            change_entropy_modified=1.5543122344752192e-15,
                            change_invariant_cubic=-5.139088546002313e-5)

        @test_allocations(semi, sol, allocs=5_000)
    end

    @trixi_testset "bbm_1d_fourier" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_fourier.jl"),
                            tspan=(0.0, 100.0),
                            l2=[6.01265210291245e-7],
                            linf=[2.3226949846799627e-7],
                            cons_error=[1.1546319456101628e-14],
                            change_waterheight=-1.1546319456101628e-14,
                            change_entropy_modified=-4.4573640600731324e-7,
                            change_invariant_cubic=-3.1871549026618595e-6)

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
