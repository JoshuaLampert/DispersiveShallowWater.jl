# DispersiveShallowWater.jl

[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JoshuaLampert.github.io/DispersiveShallowWater.jl/stable)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JoshuaLampert.github.io/DispersiveShallowWater.jl/dev/)
[![Build Status](https://github.com/JoshuaLampert/DispersiveShallowWater.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JoshuaLampert/DispersiveShallowWater.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/JoshuaLampert/DispersiveShallowWater.jl/graph/badge.svg)](https://codecov.io/gh/JoshuaLampert/DispersiveShallowWater.jl)
[![Coveralls](https://coveralls.io/repos/github/JoshuaLampert/DispersiveShallowWater.jl/badge.svg?branch=main)](https://coveralls.io/github/JoshuaLampert/DispersiveShallowWater.jl?branch=main)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/635090135.svg)](https://zenodo.org/doi/10.5281/zenodo.10034636)

[**DispersiveShallowWater.jl**](https://github.com/JoshuaLampert/DispersiveShallowWater.jl) is a [Julia](https://julialang.org/) package that implements structure-preserving numerical methods for dispersive shallow water models.
To date, it provides provably conservative, entropy-conserving and well-balanced numerical schemes for some dispersive shallow water models:

* the [BBM-BBM equations with varying bottom topography](https://iopscience.iop.org/article/10.1088/1361-6544/ac3c29),
* the [dispersive shallow water model proposed by Magnus Svärd and Henrik Kalisch](https://arxiv.org/abs/2302.09924),
* the [Serre-Green-Naghdi equations](https://arxiv.org/abs/2408.02665).

The semidiscretizations are based on summation by parts (SBP) operators, which are implemented in [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl/).
In order to obtain fully discrete schemes, the time integration methods from [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) are used to solve the resulting ordinary differential equations.
Fully discrete entropy-conservative methods can be obtained by using the [relaxation method](https://epubs.siam.org/doi/10.1137/19M1263662) provided by DispersiveShallowWater.jl.

## Installation

If you have not yet installed Julia, then you first need to [download Julia](https://julialang.org/downloads/). Please [follow the instructions for your operating system](https://julialang.org/downloads/platform/).
DispersiveShallowWater.jl works with Julia v1.9 and newer. DispersiveShallowWater.jl is a registered Julia package. Therefore, you can install it by executing the following commands from the Julia REPL
```julia
julia> using Pkg

julia> Pkg.add(["DispersiveShallowWater", "OrdinaryDiffEq", "Plots"])
```
In addition, this installs the packages [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) used for time-integration and [Plots.jl](https://github.com/JuliaPlots/Plots.jl) to visualize the results.
If you want to use other SBP operators than the default operators that DispersiveShallowWater.jl uses, then you also need [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl),
which can be installed running
```julia
julia> Pkg.add("SummationByPartsOperators")
```

## Usage

In the Julia REPL, first load the package DispersiveShallowWater.jl
```julia
julia> using DispersiveShallowWater
```

You can run a basic simulation that solves the BBM-BBM equations by executing
```julia
julia> include(default_example());
```
The result can be visualized by using the package Plots.jl
```julia
julia> using Plots
julia> plot(semi => sol)
```
The command `plot` expects a `Pair` consisting of a [`Semidiscretization`](@ref) and an `ODESolution`. The visualization can also be customized, see the [documentation](@ref visualize_results)
for more details. Other examples can be found in the subdirectory [examples/](https://github.com/JoshuaLampert/DispersiveShallowWater.jl/tree/main/examples).
A list of all examples is returned by running [`get_examples()`](@ref). You can pass the filename of one of the examples or your own simulation file to `include` in order to run it,
e.g., `include(joinpath(examples_dir(), "svaerd_kalisch_1d", "svaerd_kalisch_1d_dingemans_relaxation.jl"))`.

## Referencing

You can directly refer to DispersiveShallowWater.jl as
```bibtex
@misc{lampert2023dispersive,
  title={{D}ispersive{S}hallow{W}ater.jl: {S}tructure-preserving numerical
         methods for dispersive shallow water models},
  author={Lampert, Joshua},
  year={2023},
  month={10},
  howpublished={\url{https://github.com/JoshuaLampert/DispersiveShallowWater.jl}},
  doi={10.5281/zenodo.10034636}
}
```

## Authors

The package is developed and maintained by Joshua Lampert (University of Hamburg).
Some parts of this repository are based on parts of [Dispersive-wave-schemes-notebooks. A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations](https://github.com/ranocha/Dispersive-wave-schemes-notebooks)
by Hendrik Ranocha, Dimitrios Mitsotakis and David Ketcheson. The code structure is inspired by [Trixi.jl](https://github.com/trixi-framework/Trixi.jl/).

## License and contributing

DispersiveShallowWater.jl is published under the MIT license (see [License](https://github.com/JoshuaLampert/DispersiveShallowWater.jl/blob/main/LICENSE)). We are pleased to accept contributions from everyone, preferably in the form of a PR.
