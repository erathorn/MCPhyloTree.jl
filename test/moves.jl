@testset "tree_length_stable" begin
    tree = parsing_newick_string("(A,B,(C,(D,E)F)G)H;")
    set_binary!(tree)
    
    MCPhyloTree.number_nodes!(tree)
    
    l = tree_length(tree)
    
    swing!(tree)
    l1 = tree_length(tree)
    
    @test l == l1

end