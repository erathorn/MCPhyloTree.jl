@testset "search" begin
    tree = MCPhyloTree.parsing_newick_string("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.number_nodes!(tree)
    MCPhyloTree.set_binary!(tree)
    result = MCPhyloTree.find_num(tree,1)
    @test result.name == "A"
    leaf = MCPhyloTree.find_by_name(tree,"E")
    mid = MCPhyloTree.find_by_name(tree,"C")
    @test MCPhyloTree.find_root(leaf).name == "G"
    @test MCPhyloTree.find_root(mid).name == "G"

end
