include("../src/MCPhyloTree.jl")
using Documenter

makedocs(build   = "build",
    clean   = true,
    doctest = true,
    modules = Module[MCPhyloTree],
    sitename="MCPhyloTree",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Index" => "index.md",
        "Basics" => "basics.md",
        "Building" => "building.md",
        "Distance" => "distance.md",
        "Moves" => "moves.md",
        "Plotting" => "plotting.md"
        ]
   )

deploydocs(
     repo = "github.com/github.com/erathorn/MCPhyloTree.jl.git",
     target = "build"
     )
