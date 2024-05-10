using Documenter
using DispersiveShallowWater
using TrixiBase

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(DispersiveShallowWater, :DocTestSetup, :(using DispersiveShallowWater);
                    recursive = true)
DocMeta.setdocmeta!(TrixiBase, :DocTestSetup, :(using TrixiBase); recursive = true)

makedocs(;
         modules = [DispersiveShallowWater, TrixiBase],
         authors = "Joshua Lampert <joshua.lampert@uni-hamburg.de>",
         repo = Remotes.GitHub("JoshuaLampert", "DispersiveShallowWater.jl"),
         sitename = "DispersiveShallowWater.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://JoshuaLampert.github.io/DispersiveShallowWater.jl/stable",
                                  edit_link = "main",
                                  assets = String[],
                                  size_threshold = 1200 * 1024, # the generated .gif files can be too large
                                  size_threshold_warn = 1000 * 1024),
         pages = [
             "Home" => "index.md",
             "Overview" => "overview.md",
             "Development" => "development.md",
             "Reference" => [
                 "TrixiBase" => "ref-trixibase.md",
                 "DispersiveShallowWater" => "ref.md",
             ],
             "License" => "license.md",
         ])

deploydocs(;
           repo = "github.com/JoshuaLampert/DispersiveShallowWater.jl",
           devbranch = "main")
