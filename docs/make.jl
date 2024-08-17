using Documenter
using DispersiveShallowWater
using TrixiBase

# Dynamically replace all files in subdirectories of the source directory to include all files in these subdirectories
# This way they don't need to be listed explicitly
EQUATIONS_FILES_TO_BE_INSERTED = joinpath.(Ref("equations"),
                                           readdir(joinpath(dirname(@__DIR__), "src", "equations")))
CALLBACKS_STEP_FILES_TO_BE_INSERTED = joinpath.(Ref("callbacks_step"),
                                                readdir(joinpath(dirname(@__DIR__), "src", "callbacks_step")))

ref_path = joinpath(@__DIR__, "src", "ref.md")
lines = readlines(ref_path)
open(ref_path, "w") do io
    for line in lines
        if contains(line, "EQUATIONS_FILES_TO_BE_INSERTED")
            line = replace(line, "EQUATIONS_FILES_TO_BE_INSERTED" => ALL_FILES)
        end
        if contains(line, "CALLBACKS_STEP_FILES_TO_BE_INSERTED")
            line = replace(line, "CALLBACKS_STEP_FILES_TO_BE_INSERTED" => ALL_FILES)
        end
        println(io, line)
    end
end

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
           devbranch = "main",
           push_preview = true)
