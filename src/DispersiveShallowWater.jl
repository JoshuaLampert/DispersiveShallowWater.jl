module DispersiveShallowWater

using DiffEqBase: DiffEqBase, SciMLBase, terminate!
using Interpolations: Interpolations, linear_interpolation
using LinearAlgebra: mul!, I, Diagonal
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
                                 periodic_derivative_operator,
                                 derivative_order, integrate
import SummationByPartsOperators: grid, xmin, xmax
using TimerOutputs: TimerOutputs, TimerOutput, @timeit, print_timer, reset_timer!
@reexport using TrixiBase: TrixiBase, trixi_include

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

export boundary_condition_periodic

export initial_condition_convergence_test,
       initial_condition_manufactured, source_terms_manufactured,
       initial_condition_dingemans

export AnalysisCallback, RelaxationCallback, SummaryCallback
export tstops, errors, integrals

end
