using Documenter, MCPhyloTree

makedocs(build   = "build",
    clean   = true,
    doctest = true,
    modules = Module[MCPhyloTree],
    sitename="MCPhyloTree",
    format = Documenter.HTML(prettyurls = true),
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
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#.#"],
    target = "build"
     )
