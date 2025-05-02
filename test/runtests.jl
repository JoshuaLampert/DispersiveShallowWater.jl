using TestItems
using TestItemRunner

@run_package_tests

@testsnippet Setup begin
    include("test_util.jl")
end

@testsnippet AdditionalImports begin
    using SummationByPartsOperators: PeriodicDerivativeOperator,
                                     UniformPeriodicCoupledOperator,
                                     PeriodicUpwindOperators
    using SummationByPartsOperators: derivative_order, periodic_derivative_operator,
                                     legendre_derivative_operator,
                                     UniformPeriodicMesh1D, couple_discontinuously,
                                     upwind_operators, couple_continuously,
                                     legendre_second_derivative_operator,
                                     Mattsson2017

    using SparseArrays: sparse, SparseMatrixCSC
    using OrdinaryDiffEqTsit5: solve
    using ForwardDiff: ForwardDiff
end
