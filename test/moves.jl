
@testset "swing" begin
    tree = parsing_newick_string("((I,J)A,(K,(M,N)L)B,(C,((O,P)D,E)F)G)H;")
    set_binary!(tree)

    MCPhyloTree.number_nodes!(tree)

    l = tree_length(tree)

    swing!(tree)
    t2 = swing(tree)
    l1 = tree_length(tree)
    l2 = tree_length(t2)
    # take care of numerical instabilities
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
    # take care of numerical instabilities
    @test l == round(l2)
    @test l == round(l3)
end

@testset "NNI" begin
    tree = parsing_newick_string("((I,J)A,(K,(M,N)L)B,(C,(O,P)E)G)H;")
    set_binary!(tree)

    MCPhyloTree.number_nodes!(tree)
    tr2 = deepcopy(tree)

    res = NNI!(tr2, 10, true)

    l = tree_length(tree)
    l2 = tree_length(tr2)

    @test l == l2
    @test MCPhyloTree.RF_int(tree, tr2) == 2
end

@testset "move" begin
    tree = MCPhyloTree.parsing_newick_string("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    B = MCPhyloTree.find_by_name(tree,"B")
    A = MCPhyloTree.find_by_name(tree,"A")
    MCPhyloTree.move!(B,A,0.9)
    @test B.inc_length == 9.0
    @test A.inc_length == 1.0
end

@testset "change_edge_length" begin
    tree = MCPhyloTree.parsing_newick_string("((B:1,A:1)C:1,(D:1,E:1)F:1)G:1;")
    MCPhyloTree.number_nodes!(tree)
    MCPhyloTree.change_edge_length!(tree)
    po_list = [x.inc_length for x in MCPhyloTree.post_order(tree)]
    @test sum(po_list) != 7.0
end

@testset "reroot" begin
    tree = MCPhyloTree.parsing_newick_string("((B:1,A:1)C:1,(D:1,E:1)F:1)G:1;")
    MCPhyloTree.number_nodes!(tree)
    newtree = MCPhyloTree.reroot(tree,"C")
    @test ==MCPhyloTree.newick(newtree) == "(B:1.0,A:1.0,((D:1.0,E:1.0)F:1.0)G:1.0)C:1.0;"
end
