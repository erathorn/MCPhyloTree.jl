using Test
using MCPhyloTree
@testset "MCPhyloTree" begin
    include("basics.jl")
    println("done basics")

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
end