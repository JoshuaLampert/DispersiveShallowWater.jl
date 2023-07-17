module DispersiveShallowWater

using DiffEqBase
using LinearAlgebra: mul!, ldiv!, factorize, I, Diagonal
using PolynomialBases
using RecipesBase
using Reexport: @reexport
using Roots: AlefeldPotraShi, find_zero

using SciMLBase: CallbackSet, DiscreteCallback, ODEProblem, ODESolution
import SciMLBase: get_tmp_cache, u_modified!

@reexport using StaticArrays: SVector
using SimpleUnPack: @unpack
using SparseArrays: sparse
using SummationByPartsOperators: AbstractDerivativeOperator,
                                 UniformPeriodicMesh1D,
                                 periodic_derivative_operator,
                                 legendre_derivative_operator,
                                 couple_discontinuously,
                                 derivative_order, integrate
import SummationByPartsOperators: grid, xmin, xmax

include("util.jl")
include("boundary_conditions.jl")
include("mesh.jl")
include("solver.jl")
include("semidiscretization.jl")
include("equations/equations.jl")
include("callbacks_step/callbacks_step.jl")
include("visualization.jl")

export examples_dir, trixi_include

export BBMBBMEquations1D, BBMBBMVariableEquations1D

export waterheight_total, waterheight, velocity, momentum, energy_total, entropy,
       lake_at_rest_error

export Mesh1D, xmin, xmax, nnodes

export Solver, UpwindSolver

export Semidiscretization, semidiscretize, grid

export boundary_condition_periodic

export initial_condition_convergence_test, initial_condition_sin_bathymetry

export AnalysisCallback, RelaxationCallback
export tstops, errors, integrals

end
