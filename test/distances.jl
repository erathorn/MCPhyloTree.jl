
@testset "dist RF" begin
    tree = parsing_newick_string("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.set_binary!(tree)
    MCPhyloTree.number_nodes!(tree)
    tree2 = parsing_newick_string("((E:5.0,(B:5.0,A:5.0)C:9.0)F:5.0,D:5.0)G:5.0;")
    MCPhyloTree.set_binary!(tree)
    MCPhyloTree.number_nodes!(tree)
    @test RF(tree, tree2) == 2
end


@testset "dist RF_weighted" begin
    tree = parsing_newick_string("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.set_binary!(tree)
    MCPhyloTree.number_nodes!(tree)
    tree2 = parsing_newick_string("((E:5.0,(B:5.0,A:5.0)C:9.0)F:5.0,D:5.0)G:5.0;")
    MCPhyloTree.set_binary!(tree)
    MCPhyloTree.number_nodes!(tree)
    @test RF_weighted(tree, tree2) == 2/12
end
