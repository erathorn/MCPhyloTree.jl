@testset "upgma" begin
    D = [0. 17. 21. 31. 23. ; 17. 0. 30. 34. 21.; 21. 30. 0. 28. 39.; 31. 34. 28. 0. 43.; 23. 21. 39. 43. 0.]
    rt = upgma(D, ["a","b","c","d","e"])
    rt2 = upgma(D)
    gst = parsing_newick_string("(((a,b)u,e)v,(c,d)w)r;")
    MCPhyloTree.number_nodes!(gst)
    set_binary!(gst)
    @test RF(rt, gst) == 0
    @test RF(rt2, gst) == 0
    @test tree_height(rt) == 16.5
    @test tree_height(rt2) == 16.5
end