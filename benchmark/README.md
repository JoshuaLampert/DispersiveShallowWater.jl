# Benchmarks

This directory contains some benchmark setups using [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl).
The file `benchmarks.jl` is run by a GitHub Action leveraging [AirspeedVelocity.jl](https://github.com/MilesCranmer/AirspeedVelocity.jl)
to generate a report on the performance of the package for each pull request. If you want to run the benchmarks locally, you can do so by running
`run_benchmarks.jl`, which returns a summary `results` of the benchmark results.
