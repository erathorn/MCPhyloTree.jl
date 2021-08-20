using Documenter, MCPhyloTree

makedocs(
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
    repo = "github.com/erathorn/MCPhyloTree.jl",
    devbranch="main",
    #target = "build",
    #deps = nothing,
    #make = nothing,
    devurl = "main",
    versions = ["stable" => "v^", "v#.#", "dev" => devurl]
     )
