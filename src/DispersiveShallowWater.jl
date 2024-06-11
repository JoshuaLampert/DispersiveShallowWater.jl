module DispersiveShallowWater

using BandedMatrices: BandedMatrix
using DiffEqBase: DiffEqBase, SciMLBase, terminate!
using Interpolations: Interpolations, linear_interpolation
using LinearAlgebra: mul!, ldiv!, I, Diagonal, Symmetric, diag, lu, cholesky, cholesky!
using PolynomialBases: PolynomialBases
using Printf: @printf, @sprintf
using RecipesBase: RecipesBase, @recipe, @series
using Reexport: @reexport
using Roots: AlefeldPotraShi, find_zero

using SciMLBase: DiscreteCallback, ODEProblem, ODESolution
import SciMLBase: u_modified!

@reexport using StaticArrays: SVector
using SimpleUnPack: @unpack
using SparseArrays: sparse
using SummationByPartsOperators: AbstractDerivativeOperator,
                                 PeriodicDerivativeOperator, PeriodicUpwindOperators,
                                 UniformPeriodicCoupledOperator,
                                 DerivativeOperator, UpwindOperators,
                                 UniformCoupledOperator,
                                 periodic_derivative_operator,
                                 derivative_order, integrate, mass_matrix
import SummationByPartsOperators: grid, xmin, xmax
using TimerOutputs: TimerOutputs, print_timer, reset_timer!
@reexport using TrixiBase: trixi_include
using TrixiBase: TrixiBase, @trixi_timeit, timer

include("boundary_conditions.jl")
include("mesh.jl")
include("equations/equations.jl")
include("solver.jl")
include("semidiscretization.jl")
include("callbacks_step/callbacks_step.jl")
include("visualization.jl")
include("util.jl")

export examples_dir, get_examples, default_example, convergence_test

export BBMBBMEquations1D, BBMBBMVariableEquations1D, SvärdKalischEquations1D,
       SvaerdKalischEquations1D

export prim2prim, prim2cons, cons2prim, waterheight_total, waterheight,
       velocity, momentum, discharge, energy_total, entropy, lake_at_rest_error,
       energy_total_modified, entropy_modified

export Mesh1D, xmin, xmax, nnodes

export Solver

export Semidiscretization, semidiscretize, grid

export boundary_condition_periodic, boundary_condition_reflecting

export initial_condition_convergence_test,
       initial_condition_manufactured, source_terms_manufactured,
       initial_condition_manufactured_reflecting, source_terms_manufactured_reflecting,
       initial_condition_dingemans

export AnalysisCallback, RelaxationCallback, SummaryCallback
export tstops, errors, integrals

end
