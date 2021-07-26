using Documenter, MCPhyloTree

makedocs(build   = "build",
    clean   = true,
    doctest = true,
    modules = Module[MCPhyloTree],
    sitename="MCPhyloTree",
    format = Documenter.HTML(
        canonical = "https://erathorn.github.io/MCPhyloTree.jl/dev/"
        ),
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
    repo = "github.com/erathorn/MCPhyloTree.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    versions = ["stable" => "v^", "v#.#", devurl => "dev"]
     )
