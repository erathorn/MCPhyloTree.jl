@testset "spr" begin
    binary_tree = MCPhyloTree.parsing_newick_string("((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700,seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382);")
    MCPhyloTree.number_nodes!(binary_tree)
    MCPhyloTree.set_binary!(binary_tree)

    error_tree = MCPhyloTree.parsing_newick_string("(raccoon:19.19959):0.84600;")
    error_tree_not_binary = MCPhyloTree.parsing_newick_string("((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700,seal:12.00300,lizard:5.03):7.52973,((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382);")
    tree_binary_two = MCPhyloTree.parsing_newick_string("((raccoon:19.19959,bear:6.80041)50:0.84600,((sea_lion:11.99700,seal:12.00300)100:7.52973,((monkey:100.85930,cat:47.14069)80:20.59201,weasel:18.87953)75:2.09460)50:3.87382,dog:25.46154);")

    @test false == MCPhyloTree.check_binary(error_tree)
    @test false == MCPhyloTree.check_binary(error_tree_not_binary)
    @test_throws ErrorException MCPhyloTree.SPR(error_tree)
    @test true == MCPhyloTree.check_binary(tree_binary_two)
    spr_binary = MCPhyloTree.SPR(binary_tree)

    @test length(Set([n.num for n in post_order(spr_binary)])) == length(Set([n.num for n in post_order(binary_tree)]))

    @test round(MCPhyloTree.tree_length(binary_tree);digits=3) == round(MCPhyloTree.tree_length(spr_binary);digits=3)

    @test length(MCPhyloTree.post_order(spr_binary)) == length(MCPhyloTree.post_order(binary_tree))

    MCPhyloTree.number_nodes!(tree_binary_two)

    spr_binary_two = MCPhyloTree.SPR(tree_binary_two)

    @test length(Set([n.num for n in post_order(spr_binary_two)])) == length(Set([n.num for n in post_order(tree_binary_two)]))

    @test round(MCPhyloTree.tree_length(tree_binary_two);digits=3) == round(MCPhyloTree.tree_length(spr_binary_two);digits=3)

    @test length(MCPhyloTree.post_order(spr_binary_two)) == length(MCPhyloTree.post_order(tree_binary_two))



    spr_binary_risky = MCPhyloTree.risky_SPR(binary_tree)

    @test length(Set([n.num for n in post_order(spr_binary_risky)])) == length(Set([n.num for n in post_order(binary_tree)]))
    @test round(MCPhyloTree.tree_length(spr_binary_risky);digits=3) == round(MCPhyloTree.tree_length(spr_binary);digits=3)
    @test length(MCPhyloTree.post_order(spr_binary_risky)) == length(MCPhyloTree.post_order(binary_tree))
end