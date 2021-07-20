using Test
using MCPhyloTree
@testset "MCPhyloTree" begin
    include("basics.jl")
    println("done basics")

    include("search.jl")
    println("done search")

    include("traversal.jl")
    println("done traversal")

    include("consensus.jl")
    println("done consensus")

    include("ladderize.jl")
    println("done ladderize")

    include("pruning.jl")
    println("done pruning")

    include("spr.jl")
    println("done spr")

    include("moves.jl")
    println("done moves")

    include("building.jl")
    println("done building")
end
