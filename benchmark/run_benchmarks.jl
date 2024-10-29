using Pkg
Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
Pkg.activate(@__DIR__)
Pkg.instantiate()

include("benchmarks.jl")
tune!(SUITE)
result = run(SUITE)
