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
