using DispersiveShallowWater
using Documenter

DocMeta.setdocmeta!(DispersiveShallowWater, :DocTestSetup, :(using DispersiveShallowWater); recursive=true)

makedocs(;
    modules=[DispersiveShallowWater],
    authors="Joshua Lampert <joshua.lampert@tuhh.de>",
    repo="https://github.com/JoshuaLampert/DispersiveShallowWater.jl/blob/{commit}{path}#{line}",
    sitename="DispersiveShallowWater.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JoshuaLampert/DispersiveShallowWater.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Reference" => "ref.md",
    ],
)

#deploydocs(;
#    repo="github.com/JoshuaLampert/DispersiveShallowWater.jl",
#    devbranch="main",
#)
