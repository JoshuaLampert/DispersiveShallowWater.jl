# Changelog

DispersiveShallowWater.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.

## Changes in the v0.6 lifecycle

#### Added

- Add initial support for ForwardDiff.jl for `HyperbolicSerreGreenNaghdiEquations1D` and
  `DispersiveShallowWater.jacobian` ([#185]).

## Changes when updating to v0.6 from v0.5.x

#### Changed

- The keyword argument and function `gravity_constant` have been changed to `gravity` ([#174]).

## Changes in the v0.5 lifecycle

#### Added

- Add `LinearDispersionRelation` and documentation about dispersion ([#168]).
- Reflecting boundary conditions are added for the Svärd-Kalisch equations with `alpha = gamma = 0` ([#166]).
- Fix a bug in the upwind discretization of the `SvärdKalischEquations1D`.
- Use OrdinaryDiffEqTsit5.jl and OrdinaryDiffEqLowStorageRK.jl instead of OrdinaryDiffEq.jl in all examples to
  reduce latency ([#163]).
- Allow Fourier and periodic rational derivative operators for `BBMBBMEquations1D` and `SvärdKalischEquations1D` ([#154]).
- Add `BBMEquation1D` ([#150]).

## Changes when updating to v0.5 from v0.4.x

#### Changed

- The `BBMBBMVariableEquations1D` were removed and `BBMBBMEquations1D` now supports a `bathymetry_type` to
  choose between a flat and a variable bathymetry ([#147]).
- The default of `bathymetry_type` for the `SerreGreenNaghdiEquations1D` changed from `bathymetry_flat` to
  `bathymetry_variable` ([#147]).
- `bathymetry_type` is now a keyword argument for all equations instead of a positional argument ([#147]).
- The `initial_condition_dingemans` for the `SerreGreenNaghdiEquations1D` and `HyperbolicSerreGreenNaghdiEquations1D`
  was changed a bit to be more consistent with the other equations ([#147]).

## Changes in the v0.4 lifecycle

#### Added

- The `SerreGreenNaghdiEquations1D` were added for different types of bathymetry ([#127], [#135]).
- The `HyperbolicSerreGreenNaghdiEquations1D` were added for different types of bathymetry ([#139]).
- The abstract interface `AbstractShallowWaterEquations` was added to unify several
  systems such as the `SerreGreenNaghdiEquations1D`, the `BBMBBMEquations1D`, and the
  `SvärdKalischEquations1D` ([#127]).
- A new conversion function `prim2phys` was introduced, defaulting to `prim2prim`. `prim2phys` is the default conversion function for plotting.

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
