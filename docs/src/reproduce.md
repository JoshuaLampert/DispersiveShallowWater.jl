# How to reproduce the figures

In order to reproduce all figures used in the master thesis "Structure-preserving Numerical Methods for Dispersive Shallow Water Models" (2023) by Joshua Lampert execute the file located at the path [`DispersiveShallowWater.path_create_figures()`](@ref). From the Julia REPL, this can be done by:

```julia
julia> using DispersiveShallowWater
julia> include(DispersiveShallowWater.path_create_figures())
```
Executing this script may take a while. It will generate a folder `out/` with certain subfolders containing the figures. If you want to modify the plots or only produce a subset of plots, you can download the script [`create_figures.jl`](https://github.com/JoshuaLampert/DispersiveShallowWater.jl/blob/main/create_figures.jl), modify accordingly and run it by:

```julia
julia> include("create_figures.jl")
```
