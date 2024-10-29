using Pkg
Pkg.develop(PackageSpec(path = dirname(@__DIR__)))

include("benchmarks.jl")
tune!(SUITE)
result = run(SUITE)
