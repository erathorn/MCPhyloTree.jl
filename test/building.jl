@testset "upgma" begin
    D = [0. 17. 21. 31. 23. ; 17. 0. 30. 34. 21.; 21. 30. 0. 28. 39.; 31. 34. 28. 0. 43.; 23. 21. 39. 43. 0.]
    rt = upgma(D, ["a","b","c","d","e"])
    rt2 = upgma(D)
    gst = parsing_newick_string("(((a,b)u,e)v,(c,d)w)r;")
    MCPhyloTree.number_nodes!(gst)
    set_binary!(gst)

    @test RF(rt, gst) == 0
    @test tree_height(rt) == 16.5
    @test tree_height(rt2) == 16.5
end

@testset "NJ" begin
    D = [0. 3. 14. 12.; 3 0 13 11; 14. 13. 0. 4; 12. 11. 4. 0]
    rt = neighbor_joining(D,["a","b","d","e"])
    gst = parsing_newick_string("((a,b)c,d,e);")
    MCPhyloTree.number_nodes!(gst)
    set_binary!(gst)
    @test RF(rt, gst) == 0
    @test tree_length(rt) == 16
end

@testset "newick" begin
    tree = MCPhyloTree.parsing_newick_string("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    @test MCPhyloTree.newick(tree) == "((B:5.0,A:5.0)C:9.0,(D:5.0,E:5.0)F:5.0)G:5.0;"
end

@testset "to_covariance" begin
    tree = MCPhyloTree.parsing_newick_string("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    @test MCPhyloTree.to_covariance(tree) == [14.0 9.0 0.0 0.0; 9.0 14.0 0.0 0.0; 0.0 0.0 10.0 5.0; 0.0 0.0 5.0 10.0]
end

@testset "ladderize" begin
    tree = MCPhyloTree.parsing_newick_string("((B:5,A:5,X,Y,Z)C:9,(D:5,E:5)F:5)G:5;")
    ladder = MCPhyloTree.ladderize_tree(tree,true)
    @test MCPhyloTree.newick(ladder) == "((D:5.0,E:5.0)F:5.0,(B:5.0,A:5.0,X:1.0,Y:1.0,Z:1.0)C:9.0)G:5.0;"
