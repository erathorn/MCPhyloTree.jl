using Documenter, MCPhyloTree

makedocs(build   = "build",
    clean   = true,
    doctest = true,
    modules = Module[MCPhyloTree],
    sitename="MCPhyloTree",
    format = Documenter.HTML(),
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
    versions = ["stable" => "v^", "v#.#.#"],
    devurl = "test",
    target = "build"
     )
