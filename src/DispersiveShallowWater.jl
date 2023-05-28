module DispersiveShallowWater

using LinearAlgebra: mul!, ldiv!, factorize, I
using Reexport: @reexport

using SciMLBase: ODEProblem

@reexport using StaticArrays: SVector
using SimpleUnPack: @unpack
using SparseArrays: sparse
using SummationByPartsOperators: AbstractDerivativeOperator, periodic_derivative_operator, derivative_order
import SummationByPartsOperators: grid, xmin, xmax

include("util.jl")
include("boundary_conditions.jl")
include("mesh.jl")
include("solver.jl")
include("semidiscretization.jl")
include("equations/equations.jl")

export examples_dir, trixi_include

export BBMBBMEquations1D

export Mesh1D, xmin, xmax, nnodes

export Solver

export Semidiscretization, semidiscretize, grid

export boundary_condition_periodic

export initial_condition_convergence_test

end
