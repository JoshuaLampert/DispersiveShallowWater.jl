module TestUnit

using Test
using DispersiveShallowWater
using SummationByPartsOperators: PeriodicDerivativeOperator, UniformPeriodicCoupledOperator
using SummationByPartsOperators: derivative_order, legendre_derivative_operator,
                                 UniformPeriodicMesh1D, couple_discontinuously
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
    mesh = Mesh1D(-1, 1, 10)
    p = 3
    solver = @test_nowarn Solver(mesh, p)
    @test solver.D isa PeriodicDerivativeOperator
    @test solver.D2 isa PeriodicDerivativeOperator
    @test derivative_order(solver.D) == 1
    @test derivative_order(solver.D2) == 2
    @test grid(solver) == grid(solver.D) == grid(solver.D2)
    @test real(solver) == Float64
    @test_nowarn show(stdout, solver)

    Dop = legendre_derivative_operator(-1.0, 1.0, p + 1)
    sbp_mesh = UniformPeriodicMesh1D(-1.0, 1.0, 512 ÷ (p + 1))
    D = couple_discontinuously(Dop, sbp_mesh)
    D_pl = couple_discontinuously(Dop, sbp_mesh, Val(:plus))
    D_min = couple_discontinuously(Dop, sbp_mesh, Val(:minus))
    D2 = sparse(D_pl) * sparse(D_min)
    solver = @test_nowarn Solver(D, D2)
    @test solver.D isa UniformPeriodicCoupledOperator
    @test solver.D2 isa SparseMatrixCSC
  end
end

end # module
