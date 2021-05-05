@testset "swing" begin
    tree = parsing_newick_string("((I,J)A,(K,(M,N)L)B,(C,((O,P)D,E)F)G)H;")
    set_binary!(tree)
    
    MCPhyloTree.number_nodes!(tree)
    
    l = tree_length(tree)
    
    swing!(tree)
    t2 = swing(tree)
    l1 = tree_length(tree)
    l2 = tree_length(t2)
    @test l == l1
    @test l == l2
end

@testset "slide" begin
    tree = parsing_newick_string("((I,J)A,(K,(M,N)L)B,(C,((O,P)D,E)F)G)H;")
    set_binary!(tree)
    
    MCPhyloTree.number_nodes!(tree)
    
    l = tree_length(tree)
    
    slide!(tree)
    l2 = tree_length(tree)

    t2 = slide(tree)
    l3 = tree_length(t2)

    @test l == l2
    @test l == l3
end