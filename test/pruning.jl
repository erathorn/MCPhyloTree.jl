
tree = MCPhyloTree.parsing_newick_string("(A,B,(C,(D,E)F)G)H;")
MCPhyloTree.number_nodes!(tree)

@testset "Pruning leaf" begin
    node_list = ["C"]
    prune_tree = MCPhyloTree.prune_tree(tree, node_list)
    @test newick(prune_tree) == newick(MCPhyloTree.parsing_newick_string("(A,B,((D,E)F)G)H;"))
end

@testset "Pruning root" begin
    node_list = ["H"]
    @test_throws ArgumentError MCPhyloTree.prune_tree(tree, node_list)
end

@testset "Pruning misc" begin
    node_list1 = ["B", "C"]
    node_list2 = ["F"]
    node_list3 = ["A", "C", "G"]

    prune_tree1 = MCPhyloTree.prune_tree(tree, node_list1)
    prune_tree2 = MCPhyloTree.prune_tree(tree, node_list2)
    prune_tree3 = MCPhyloTree.prune_tree(tree, node_list3)

    @test newick(prune_tree1) == newick(MCPhyloTree.parsing_newick_string("(A,((D,E)F)G)H;"))
    @test newick(prune_tree2) == newick(MCPhyloTree.parsing_newick_string("(A,B,(C)G)H;"))
    @test newick(prune_tree3) == newick(MCPhyloTree.parsing_newick_string("(B)H;"))
end

@testset "Pruning inplace" begin
    node_list1 = ["B", "C"]
    MCPhyloTree.prune_tree!(tree, node_list1)
    @test newick(tree) ==  newick(MCPhyloTree.parsing_newick_string("(A,((D,E)F)G)H;"))
end
