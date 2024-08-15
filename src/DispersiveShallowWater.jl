"""
    DispersiveShallowWater

**DispersiveShallowWater.jl** is a Julia package that implements structure-preserving numerical methods for dispersive shallow water models.
It provides provably conservative, entropy-conserving and well-balanced numerical schemes for some dispersive shallow water models.

The semidiscretizations are based on summation-by-parts (SBP) operators, which are implemented in SummationByPartsOperators.jl. To
obtain fully discrete schemes, the time integration methods from OrdinaryDiffEq.jl are used to solve the resulting ordinary differential equations.
Fully discrete entropy-conservative methods can be obtained by using the relaxation method provided by DispersiveShallowWater.jl.

See also: [DispersiveShallowWater.jl](https://github.com/JoshuaLampert/DispersiveShallowWater.jl)
"""
module DispersiveShallowWater

using BandedMatrices: BandedMatrix
using DiffEqBase: DiffEqBase, SciMLBase, terminate!
using Interpolations: Interpolations, linear_interpolation
using LinearAlgebra: mul!, ldiv!, I, Diagonal, Symmetric, diag, lu, cholesky, cholesky!,
                     issuccess
using PolynomialBases: PolynomialBases
using Printf: @printf, @sprintf
using RecipesBase: RecipesBase, @recipe, @series
using RecursiveArrayTools: ArrayPartition
using Reexport: @reexport
using Roots: AlefeldPotraShi, find_zero

using SciMLBase: DiscreteCallback, ODEProblem, ODESolution
import SciMLBase: u_modified!

@reexport using StaticArrays: SVector
using SimpleUnPack: @unpack
using SparseArrays: sparse, issparse
using SummationByPartsOperators: SummationByPartsOperators,
                                 AbstractDerivativeOperator,
                                 PeriodicDerivativeOperator, PeriodicUpwindOperators,
                                 UniformPeriodicCoupledOperator,
                                 DerivativeOperator, UpwindOperators,
                                 UniformCoupledOperator,
                                 FourierDerivativeOperator,
                                 periodic_derivative_operator,
                                 derivative_order, integrate, mass_matrix,
                                 scale_by_mass_matrix!
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

export AbstractShallowWaterEquations,
       BBMBBMEquations1D, BBMBBMVariableEquations1D,
       SvärdKalischEquations1D, SvaerdKalischEquations1D,
       SerreGreenNaghdiEquations1D

export prim2prim, prim2cons, cons2prim,
       waterheight_total, waterheight,
       velocity, momentum, discharge,
       gravity_constant,
       bathymetry, still_water_surface,
       energy_total, entropy, lake_at_rest_error,
       energy_total_modified, entropy_modified

export Mesh1D, xmin, xmax, nnodes

export Solver

export Semidiscretization, semidiscretize, grid

export boundary_condition_periodic, boundary_condition_reflecting

export bathymetry_flat

export initial_condition_convergence_test,
       initial_condition_manufactured, source_terms_manufactured,
       initial_condition_manufactured_reflecting, source_terms_manufactured_reflecting,
       initial_condition_dingemans

export AnalysisCallback, RelaxationCallback, SummaryCallback
export tstops, errors, integrals

end
