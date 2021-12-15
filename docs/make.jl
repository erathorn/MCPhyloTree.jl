using Documenter, MCPhyloTree

makedocs(
    modules = MCPhyloTree,
    sitename="MCPhyloTree.jl",
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
    devbranch="main",
    versions = ["stable" => "v^", "v#.#.#", devurl => devurl]
     )
