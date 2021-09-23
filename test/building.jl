@testset "upgma" begin
    D = [
        0.0 17.0 21.0 31.0 23.0
        17.0 0.0 30.0 34.0 21.0
        21.0 30.0 0.0 28.0 39.0
        31.0 34.0 28.0 0.0 43.0
        23.0 21.0 39.0 43.0 0.0
    ]
    rt = upgma(D, ["a", "b", "c", "d", "e"])
    rt2 = upgma(D)
    gst = ParseNewick("(((a,b)u,e)v,(c,d)w)r;")

    @test RF(rt, gst) == 0
    @test tree_height(rt) == 16.5
    @test tree_height(rt2) == 16.5
end

@testset "NJ" begin
    D = [0.0 3.0 14.0 12.0; 3 0 13 11; 14.0 13.0 0.0 4; 12.0 11.0 4.0 0]
    rt = neighbor_joining(D, ["a", "b", "d", "e"])
    gst = ParseNewick("((a,b)c,d,e);")
    @test RF(rt, gst) == 0
    @test tree_length(rt) == 16
end

@testset "newick" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    @test newick(tree) == "((B:5.0,A:5.0)C:9.0,(D:5.0,E:5.0)F:5.0)G:5.0;"
end

@testset "to_covariance" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    @test to_covariance(tree) ==
          [14.0 9.0 0.0 0.0; 9.0 14.0 0.0 0.0; 0.0 0.0 10.0 5.0; 0.0 0.0 5.0 10.0]
end

@testset "ladderize" begin
    tree = ParseNewick("((B:5,A:5,X,Y,Z)C:9,(D:5,E:5)F:5)G:5;")
    ladder = ladderize_tree(tree, true)
    @test newick(ladder) ==
          "((D:5.0,E:5.0)F:5.0,(B:5.0,A:5.0,X:1.0,Y:1.0,Z:1.0)C:9.0)G:5.0;"
end

@testset "to_df" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    target_arr = [
        0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0
        5.0 5.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 5.0 5.0 0.0 0.0
        0.0 0.0 9.0 0.0 0.0 5.0 0.0
    ]
    target_names = ["B", "A", "C", "D", "E", "F", "G"]
    res_arr, res_names = to_df(tree)
    @test all(target_arr .== res_arr)
    @test all(target_names .== res_names)
end


@testset "from_df" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    target_arr = [
        0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0
        5.0 5.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 5.0 5.0 0.0 0.0
        0.0 0.0 9.0 0.0 0.0 5.0 0.0
    ]
    target_names = ["B", "A", "C", "D", "E", "F", "G"]
    res_tree = from_df(target_arr, target_names)
    @test RF(res_tree, tree) == 0
end


@testset "leave_incidence_matrix" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    lm = leave_incidence_matrix(tree)
    target = [
        0.0 1.0 0.0 0.0 1.0 0.0
        1.0 0.0 0.0 0.0 1.0 0.0
        0.0 0.0 1.0 0.0 0.0 1.0
        0.0 0.0 0.0 1.0 0.0 1.0
    ]
    @test all(lm .== target)
end
