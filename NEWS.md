# Changelog

DispersiveShallowWater.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.


## Changes in the v0.4 lifecycle

#### Added
- The `SerreGreenNaghdiEquations1D` were added for different types of bathymetry ([#127], [#135]).
- The `HyperbolicSerreGreenNaghdiEquations1D` were added for different types of bathymetry ([#139]).
- The abstract interface `AbstractShallowWaterEquations` was added to unify several
  systems such as the `SerreGreenNaghdiEquations1D`, the `BBMBBMEquations1D`, and the
  `SvärdKalischEquations1D` ([#127]).

## Changes when updating to v0.4 from v0.3.x

#### Changed
- Use `ArrayPartition` from RecursiveArrayTools.jl to store the solution of the `ODEProblem` ([#118]).

## Changes in the v0.3 lifecycle

#### Added
- Add possibility to pass vector of `Ns` to `convergence_test` ([#113]).
- Performance improvements by using factorized matrices for linear systems solves ([#108], [#112], [#114]).
- Reflecting boundary conditions are added for the BBM-BBM equations ([#104], [#109]).
- Fix for the `BBMBBMVariableEquations1D`, where the still water surface was neglected leading
  to a bug in the Dingemans setup ([#91]).

## Changes when updating to v0.3 from v0.2.x

#### Changed
- Add keyword argument `start_from` when plotting `AnalysisCallback` ([#87]).
- Manufactured solution for Svärd-Kalisch equations uses a variable bathymetry ([#84]).

## Changes in the v0.2 lifecycle

#### Added
- Add `SummaryCallback` ([#75]).

## Changes when updating to v0.2 from v0.1.x

#### Changed
- The code from the master thesis of Joshua Lampert was separated ([#69]).
- Add support for source terms ([#65]).
- A higher order interpolation is used when plotting the solution at a value `x` outside
  the grid ([#64]).
