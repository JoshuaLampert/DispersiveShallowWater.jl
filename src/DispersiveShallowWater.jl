module DispersiveShallowWater

using DiffEqBase
using LinearAlgebra: mul!, ldiv!, factorize, I, Diagonal
using PolynomialBases
using Printf: @printf, @sprintf
using RecipesBase
using Reexport: @reexport
using Roots: AlefeldPotraShi, find_zero

using SciMLBase: CallbackSet, DiscreteCallback, ODEProblem, ODESolution
import SciMLBase: get_tmp_cache, u_modified!

@reexport using StaticArrays: SVector
using SimpleUnPack: @unpack
using SparseArrays: sparse
using SummationByPartsOperators: AbstractDerivativeOperator,
                                 PeriodicDerivativeOperator, PeriodicUpwindOperators,
                                 periodic_derivative_operator,
                                 derivative_order, integrate
import SummationByPartsOperators: grid, xmin, xmax

include("util.jl")
include("boundary_conditions.jl")
include("mesh.jl")
include("solver.jl")
include("equations/equations.jl")
include("semidiscretization.jl")
include("callbacks_step/callbacks_step.jl")
include("visualization.jl")

export examples_dir, trixi_include

export BBMBBMEquations1D, BBMBBMVariableEquations1D, SvärdKalischEquations1D,
       SvaerdKalischEquations1D

export prim2prim, prim2cons, cons2prim, waterheight_total, waterheight,
       velocity, momentum, discharge, energy_total, entropy, lake_at_rest_error,
       energy_total_modified, entropy_modified

export Mesh1D, xmin, xmax, nnodes

export Solver

export Semidiscretization, semidiscretize, grid

export boundary_condition_periodic

export initial_condition_convergence_test, initial_condition_sin_bathymetry,
       initial_condition_dingemans

export AnalysisCallback, RelaxationCallback
export tstops, errors, integrals

end
