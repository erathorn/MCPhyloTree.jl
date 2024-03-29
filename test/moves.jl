
@testset "swing" begin
    tree = ParseNewick("((I,J)A,(K,(M,N)L)B,(C,((O,P)D,E)F)G)H;")

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
    tree = ParseNewick("((I,J)A,(K,(M,N)L)B,(C,((O,P)D,E)F)G)H;")

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
    tree = ParseNewick("((I,J)A,(K,(M,N)L)B,(C,(O,P)E)G)H;")

    tr2 = deepcopy(tree)

    res = NNI!(tr2, 10, true)

    l = tree_length(tree)
    l2 = tree_length(tr2)

    @test l == l2
    @test MCPhyloTree.RF_int(tree, tr2)[1] == 2

    ## NNI with root
    @test NNI!(tr2, tr2.num, true) == 0

    ## NNI without target direction
    tr2 = deepcopy(tree)

    res = NNI!(tr2, 10)

    l = tree_length(tree)
    l2 = tree_length(tr2)

    @test l == l2
    @test MCPhyloTree.RF_int(tree, tr2)[1] == 2
    

end

@testset "move" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    B = find_by_name(tree,"B")
    A = find_by_name(tree,"A")
    MCPhyloTree.move!(B,A,0.9)
    @test B.inc_length == 9.0
    @test A.inc_length == 1.0
end

@testset "change_edge_length" begin
    tree = ParseNewick("((B:1,A:1)C:1,(D:1,E:1)F:1)G:1;")
    change_edge_length!(tree)
    po_list = [x.inc_length for x in post_order(tree)]
    @test sum(po_list) != 7.0
end

@testset "reroot" begin
    tree = ParseNewick("((B:1,A:1)C:1,(D:1,E:1)F:1)G:1;")
    newtree = MCPhyloTree.reroot(tree,"C")
    @test newick(newtree) == "(B:1.0,A:1.0,((D:1.0,E:1.0)F:1.0)G:1.0)C:1.0;"
end

@testset "SPR" begin
    tree2 = ParseNewick("(A,B,(C,(D,E)F)G)H;")
    subtree = find_by_name(tree2,"F")
    target = find_by_name(tree2,"A")
    nutree = MCPhyloTree.perform_spr(tree2,subtree,target)
    @test newick(nutree) == "(B:1.0,C:1.0,((D:1.0,E:1.0)F:1.0,A:1.0)G:1.0)H:1.0;"
end

@testset "NNI Matrix" begin
    target = [
        1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  1.0  1.0  1.0  1.0
        0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0  1.0  1.0
        0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0  0.0  1.0  1.0  1.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0  1.0  0.0  0.0  1.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  1.0  0.0  0.0  1.0  1.0  1.0
    ]
    lm1 = [
        1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0
        0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0  1.0
        0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0
    ]

    res = MCPhyloTree.NNI_int!(lm1, 11, 10)
    @test res == 1
    @test all(lm1 .== target)

    lm1 = [
        1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0
        0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0  1.0
        0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0
    ]
    n_lm =deepcopy(lm1)
    while NNI!(n_lm) != 1
    end
    @test RF(n_lm, lm1) == 2
end

