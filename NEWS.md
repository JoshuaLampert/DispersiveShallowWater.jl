# Changelog

DispersiveShallowWater.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.


## Changes in the v0.4 lifecycle

#### Added
- The `SerreGreenNaghdiEquations1D` were added for different types of bathymetry.
- The abstract interface `AbstractShallowWaterEquations` was added to unify several
  systems such as the `SerreGreenNaghdiEquations1D`, the `BBMBBMEquations1D`, and the
  `SvärdKalischEquations1D`.
