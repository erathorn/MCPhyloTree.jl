using Test
using MCPhyloTree
@testset "MCPhyloTree" begin
    include("basics.jl")
    include("consensus.jl")
    include("ladderize.jl")
    include("pruning.jl")
    include("spr.jl")
    include("moves.jl")
end