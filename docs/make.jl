using Documenter
using DispersiveShallowWater
using TrixiBase
using Changelog: Changelog

# Dynamically set all files in subdirectories of the source directory to include all files in these subdirectories
# This way they don't need to be listed explicitly
EQUATIONS_FILES = joinpath.(Ref("equations"),
                            readdir(joinpath(dirname(@__DIR__), "src",
                                             "equations")))
CALLBACKS_STEP_FILES = joinpath.(Ref("callbacks_step"),
                                 readdir(joinpath(dirname(@__DIR__), "src",
                                                  "callbacks_step")))

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(DispersiveShallowWater, :DocTestSetup, :(using DispersiveShallowWater);
                    recursive = true)
DocMeta.setdocmeta!(TrixiBase, :DocTestSetup, :(using TrixiBase); recursive = true)

# Copy some files from the repository root directory to the docs and modify them
# as necessary
# Based on: https://github.com/ranocha/SummationByPartsOperators.jl/blob/0206a74140d5c6eb9921ca5021cb7bf2da1a306d/docs/make.jl#L27-L41
open(joinpath(@__DIR__, "src", "code_of_conduct.md"), "w") do io
    # Point to source file
    println(io,
            """
            ```@meta
            EditURL = "https://github.com/JoshuaLampert/DispersiveShallowWater.jl/blob/main/CODE_OF_CONDUCT.md"
            ```
            """)
    # Write the modified contents
    println(io, "# [Code of Conduct](@id code-of-conduct)")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"))
        println(io, "> ", line)
    end
end

open(joinpath(@__DIR__, "src", "contributing.md"), "w") do io
    # Point to source file
    println(io,
            """
            ```@meta
            EditURL = "https://github.com/JoshuaLampert/DispersiveShallowWater.jl/blob/main/CONTRIBUTING.md"
            ```
            """)
    # Write the modified contents
    for line in eachline(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"))
        line = replace(line, "[LICENSE](LICENSE)" => "[License](@ref)")
        println(io, line)
    end
end

# Create changelog
Changelog.generate(Changelog.Documenter(),                           # output type
                   joinpath(@__DIR__, "..", "NEWS.md"),              # input file
                   joinpath(@__DIR__, "src", "changelog_tmp.md");    # output file
                   repo = "JoshuaLampert/DispersiveShallowWater.jl", # default repository for links
                   branch = "main",)
# Fix edit URL of changelog
open(joinpath(@__DIR__, "src", "changelog.md"), "w") do io
    for line in eachline(joinpath(@__DIR__, "src", "changelog_tmp.md"))
        if startswith(line, "EditURL")
            line = "EditURL = \"https://github.com/JoshuaLampert/DispersiveShallowWater.jl/blob/main/NEWS.md\""
        end
        println(io, line)
    end
end

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
         pages = ["Home" => "index.md",
             "Overview" => "overview.md",
             "Development" => "development.md",
             "Reference" => [
                 "TrixiBase" => "ref-trixibase.md",
                 "DispersiveShallowWater" => "ref.md",
             ],
             "Changelog" => "changelog.md",
             "Contributing" => "contributing.md",
             "Code of Conduct" => "code_of_conduct.md",
             "License" => "license.md"])

deploydocs(;
           repo = "github.com/JoshuaLampert/DispersiveShallowWater.jl",
           devbranch = "main",
           push_preview = true)
