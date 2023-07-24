module TestUnit

using Test
using DispersiveShallowWater
using SummationByPartsOperators: PeriodicDerivativeOperator, UniformPeriodicCoupledOperator,
                                 PeriodicUpwindOperators
using SummationByPartsOperators: derivative_order, periodic_derivative_operator,
                                 legendre_derivative_operator,
                                 UniformPeriodicMesh1D, couple_discontinuously,
                                 upwind_operators
using SparseArrays: sparse, SparseMatrixCSC

@testset "Unit tests" begin
    @testset "Mesh1D" begin
        mesh = @test_nowarn Mesh1D(-1, 1, 10)
        @test ndims(mesh) == 1
        @test xmin(mesh) == -1
        @test xmax(mesh) == 1
        @test nnodes(mesh) == 10
        @test real(mesh) == Int64
        @test_nowarn show(stdout, mesh)
    end

    @testset "Solver" begin
        mesh = Mesh1D(-1.0, 1.0, 10)
        p = 3
        solver = @test_nowarn Solver(mesh, p)
        @test solver.D1 isa PeriodicDerivativeOperator
        @test solver.D2 isa PeriodicDerivativeOperator
        @test derivative_order(solver.D1) == 1
        @test derivative_order(solver.D2) == 2
        @test grid(solver) == grid(solver.D1) == grid(solver.D2)
        @test real(solver) == Float64
        @test_nowarn show(stdout, solver)

        Dop = legendre_derivative_operator(-1.0, 1.0, p + 1)
        sbp_mesh = UniformPeriodicMesh1D(-1.0, 1.0, 512 ÷ (p + 1))
        central = couple_discontinuously(Dop, sbp_mesh)
        minus = couple_discontinuously(Dop, sbp_mesh, Val(:minus))
        plus = couple_discontinuously(Dop, sbp_mesh, Val(:plus))
        D2 = sparse(plus) * sparse(minus)
        solver = @test_nowarn Solver(central, D2)
        @test solver.D1 isa UniformPeriodicCoupledOperator
        @test solver.D2 isa SparseMatrixCSC
        D1 = PeriodicUpwindOperators(minus, central, plus)
        solver = @test_nowarn Solver(D1, D2)
        @test solver.D1 isa PeriodicUpwindOperators
        @test solver.D2 isa SparseMatrixCSC

        D1 = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                              accuracy_order = p, xmin = -1.0, xmax = 1.0, N = 10)
        D2 = sparse(D1.plus) * sparse(D1.minus)
        solver = @test_nowarn Solver(D1, D2)
        @test solver.D1 isa PeriodicUpwindOperators
        @test solver.D2 isa SparseMatrixCSC
        @test derivative_order(solver.D1) == 1
    end

    @testset "BBMBBMEquations1D" begin
        equations = @test_nowarn BBMBBMEquations1D(gravity_constant = 9.81, D = 2.0)
        u = [42.0, 2.0]
        @test waterheight_total(u, equations) == 42.0
        @test waterheight(u, equations) == 44.0
        @test velocity(u, equations) == 2.0
        @test momentum(u, equations) == 88.0
        @test isapprox(energy_total(u, equations), 17480.84)
        @test_nowarn show(stdout, equations)
    end

    @testset "BBMBBMVariableEquations1D" begin
        equations = @test_nowarn BBMBBMVariableEquations1D(gravity_constant = 9.81)
        u = [42.0, 2.0, 2.0]
        @test waterheight_total(u, equations) == 42.0
        @test waterheight(u, equations) == 44.0
        @test velocity(u, equations) == 2.0
        @test momentum(u, equations) == 88.0
        @test isapprox(energy_total(u, equations), 17480.84)
        @test_nowarn show(stdout, equations)
    end

    @testset "SvaerdKalischEquations1D" begin
        equations = @test_nowarn SvaerdKalischEquations1D(gravity_constant = 9.81,
                                                          alpha = 0.0004040404040404049,
                                                          beta = 0.49292929292929294,
                                                          gamma = 0.15707070707070708)
        u = [42.0, 2.0, 2.0]
        @test waterheight_total(u, equations) == 42.0
        @test waterheight(u, equations) == 44.0
        @test velocity(u, equations) == 2.0
        @test momentum(u, equations) == 88.0
        @test isapprox(energy_total(u, equations), 17480.84)
        @test_nowarn show(stdout, equations)
    end
end

end # module
