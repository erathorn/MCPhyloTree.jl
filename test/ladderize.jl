
@testset "ladderize" begin
    start_tree = MCPhyloTree.parsing_newick_string("((A,(B,(C,(D,E)))),(F,(G,H)))")
    descending_tree = MCPhyloTree.parsing_newick_string("(((((D,E),C),B),A),((G,H),F))")
    ascending_tree = MCPhyloTree.parsing_newick_string("((F,(G,H)),(A,(B,(C,(D,E)))))")

    start_newick = newick(start_tree)
    desc_newick = newick(descending_tree)
    asc_newick = newick(ascending_tree)

    ladderized_asc_tree = MCPhyloTree.ladderize_tree(start_tree)
    @test newick(ladderized_asc_tree) == asc_newick

    ladderized_desc_tree = MCPhyloTree.ladderize_tree(start_tree,false)
    @test newick(ladderized_desc_tree) == desc_newick

    MCPhyloTree.ladderize_tree!(start_tree)
    @test newick(start_tree) == asc_newick

    start_tree = MCPhyloTree.parsing_newick_string("((A,(B,(C,(D,E)))),(F,(G,H)))")
    MCPhyloTree.ladderize_tree!(start_tree, false)
    @test newick(start_tree) == desc_newick
end