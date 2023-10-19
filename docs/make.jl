using DispersiveShallowWater
using Documenter

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(DispersiveShallowWater, :DocTestSetup, :(using DispersiveShallowWater);
                    recursive = true)

makedocs(;
         modules = [DispersiveShallowWater],
         authors = "Joshua Lampert <joshua.lampert@uni-hamburg.de>",
         repo = "https://github.com/JoshuaLampert/DispersiveShallowWater.jl/blob/{commit}{path}#{line}",
         sitename = "DispersiveShallowWater.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://JoshuaLampert.github.io/DispersiveShallowWater.jl/stable",
                                  edit_link = "main",
                                  assets = String[]),
         pages = [
             "Home" => "index.md",
             "Overview" => "overview.md",
             "Reproduce figures" => "reproduce.md",
             "Reference" => "ref.md",
             "License" => "license.md",
         ])

deploydocs(;
           repo = "github.com/JoshuaLampert/DispersiveShallowWater.jl",
           devbranch = "main")
