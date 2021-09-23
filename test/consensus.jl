
@testset "find_common_clusters" begin
    ref_tree = ParseNewick("((A,(B,(C,(D,E)))),(F,(G,H)));")
    tree2 = ParseNewick("((G,(C,(A,(F,E)))),(B,(D,H)));")

    A, B, C, D, E, F, G, H = find_by_name(tree2, "A"), find_by_name(tree2, "B"),
                             find_by_name(tree2, "C"), find_by_name(tree2, "D"),
                             find_by_name(tree2, "E"), find_by_name(tree2, "F"),
                             find_by_name(tree2, "G"), find_by_name(tree2, "H")
     expected_dict = Dict{Int64,Tuple{Bool,Union{Missing, Float64}}}([(D.mother.num, (false, missing)), (C.mother.num, (false, missing)), (B.mother.num, (false, missing)),
                           (A.mother.num, (false, missing)), (G.mother.num, (false, missing)), (F.mother.num, (false, missing)),
                           (G.mother.mother.num, (true, 1.0)), (A.num, (true, 1.0)), (B.num, (true, 1.0)), (C.num, (true, 1.0)),
                           (D.num, (true, 1.0)), (E.num, (true, 1.0)), (F.num, (true, 1.0)), (G.num, (true, 1.0)), (H.num, (true, 1.0))])
    @test isequal(MCPhyloTree.find_common_clusters(ref_tree, tree2) , expected_dict)

    tree3 = ParseNewick("((A,(C,(D,(B,E)))),(G,(F,H)));")
    
    A, B, C, D, E, F, G, H = find_by_name(tree3, "A"), find_by_name(tree3, "B"),
                             find_by_name(tree3, "C"), find_by_name(tree3, "D"),
                             find_by_name(tree3, "E"), find_by_name(tree3, "F"),
                             find_by_name(tree3, "G"), find_by_name(tree3, "H")
     expected_dict = Dict([(D.mother.num, (false, missing)), (C.mother.num, (true, 1.0)), (B.mother.num, (false, missing)),
                           (A.mother.num, (true, 1.0)), (G.mother.num, (true, 1.0)), (F.mother.num, (false, missing)),
                           (A.mother.mother.num, (true, 1.0)), (A.num, (true, 1.0)), (B.num, (true, 1.0)), (C.num, (true, 1.0)),
                           (D.num, (true, 1.0)), (E.num, (true, 1.0)), (F.num, (true, 1.0)), (G.num, (true, 1.0)), (H.num, (true, 1.0))])
    @test isequal(MCPhyloTree.find_common_clusters(ref_tree, tree3), expected_dict)

    tree4 = ParseNewick("((A,(B,(C,(D,E)))),(F,(G,H)));")
    
    A, B, C, D, E, F, G, H = find_by_name(tree4, "A"), find_by_name(tree4, "B"),
                             find_by_name(tree4, "C"), find_by_name(tree4, "D"),
                             find_by_name(tree4, "E"), find_by_name(tree4, "F"),
                             find_by_name(tree4, "G"), find_by_name(tree4, "H")
    expected_dict = Dict([(D.mother.num, (true, 1.0)), (C.mother.num, (true, 1.0)), (B.mother.num, (true, 1.0)),
                          (A.mother.num, (true, 1.0)), (G.mother.num, (true, 1.0)), (F.mother.num, (true, 1.0)),
                          (A.mother.mother.num, (true, 1.0)), (A.num, (true, 1.0)), (B.num, (true, 1.0)), (C.num, (true, 1.0)),
                          (D.num, (true, 1.0)), (E.num, (true, 1.0)), (F.num, (true, 1.0)), (G.num, (true, 1.0)), (H.num, (true, 1.0))])
    @test isequal(MCPhyloTree.find_common_clusters(ref_tree, tree4), expected_dict)

    tree5 = ParseNewick("((G,(X,(A,(F,E)))),(B,(D,H)));")
    
    @test_throws ArgumentError MCPhyloTree.find_common_clusters(ref_tree, tree5)

    tree6 = ParseNewick("(X,(G,(C,(A,(F,E)))),(B,(D,H)));")
 
    @test_throws ArgumentError MCPhyloTree.find_common_clusters(ref_tree, tree6)
end

@testset "one_way_compatible" begin
    tree = ParseNewick("((A,C,E),(B,D));")
    tree2 = ParseNewick("((A,C),(B,D,E));")
    expected_tree = ParseNewick("(A,C,E,(B,D));")
  
    @test newick(MCPhyloTree.one_way_compatible(tree, tree2)) == newick(expected_tree)
end

@testset "get_leaf_ranks" begin
    ref_tree = ParseNewick("((A,(B,(C,(D,E)))),(F,(G,H)));")
    nodes = post_order(ref_tree)
    @test MCPhyloTree.get_leaf_ranks(nodes) == Dict([("A", 1), ("B", 2), ("C", 3),
                                                 ("D", 4), ("E", 5), ("F", 6),
                                                 ("G", 7), ("H", 8)])
end

@testset "order_tree!" begin
    tree = ParseNewick("(A,B,(C,(D,E)F)G)H;")
    
    A, B, C, D, E, F, G, H = find_by_name(tree, "A"), find_by_name(tree, "B"),
                             find_by_name(tree, "C"), find_by_name(tree, "D"),
                             find_by_name(tree, "E"), find_by_name(tree, "F"),
                             find_by_name(tree, "G"), find_by_name(tree, "H")
    cluster_start_indeces = Dict([(A, 3), (B, 7), (C, 2), (D, 8),
                                  (E, 5), (F, 1), (G, 4), (H, 6)])
    ordered_tree = ParseNewick("(A,((E,D)F,C)G,B)H;")
    @test MCPhyloTree.order_tree!(tree, cluster_start_indeces) == [A, E, D, C, B]
    @test newick(tree) == newick(ordered_tree)
end

@testset "max/min_leaf_rank" begin
    tree = ParseNewick("(A,B,(C,(D,E)F)G)H;")
    F, G, H, A = find_by_name(tree, "F"), find_by_name(tree, "G"),
                 find_by_name(tree, "H"), find_by_name(tree, "A")
    leaf_ranks = Dict([("A", 1), ("B", 5), ("C", 2), ("D", 4), ("E", 3)])

    @testset "min_leaf_rank" begin
        @test MCPhyloTree.min_leaf_rank(leaf_ranks, F) == 3
        @test MCPhyloTree.min_leaf_rank(leaf_ranks, G) == 2
        @test MCPhyloTree.min_leaf_rank(leaf_ranks, H) == 1
        @test MCPhyloTree.min_leaf_rank(leaf_ranks, A) == 1
    end

    @testset "max_leaf_rank" begin
        @test MCPhyloTree.max_leaf_rank(leaf_ranks, F) == 4
        @test MCPhyloTree.max_leaf_rank(leaf_ranks, G) == 4
        @test MCPhyloTree.max_leaf_rank(leaf_ranks, H) == 5
        @test MCPhyloTree.max_leaf_rank(leaf_ranks, A) == 1
    end
end

@testset "x_left_right" begin
    tree = ParseNewick("(A,B,(C,(D,E)F)G)H;")
    
    A, B, C, D, E, F, G, H = find_by_name(tree, "A"), find_by_name(tree, "B"),
                             find_by_name(tree, "C"), find_by_name(tree, "D"),
                             find_by_name(tree, "E"), find_by_name(tree, "F"),
                             find_by_name(tree, "G"), find_by_name(tree, "H")

    @testset "x_left" begin
        @test MCPhyloTree.x_left(A) == (H, [A,H])
        @test MCPhyloTree.x_left(B) == (B, [B,H])
        @test MCPhyloTree.x_left(C) == (G, [C,G,H])
        @test MCPhyloTree.x_left(D) == (F, [D,F,G])
        @test MCPhyloTree.x_left(E) == (E, [E,F])
    end

    @testset "x_right" begin
        @test MCPhyloTree.x_right(A) == (A, [A,H])
        @test MCPhyloTree.x_right(B) == (B, [B,H])
        @test MCPhyloTree.x_right(C) == (C, [C,G])
        @test MCPhyloTree.x_right(D) == (D, [D,F])
        @test MCPhyloTree.x_right(E) == (H, [E,F,G,H])
    end
end

@testset "majority_consensus_tree" begin
    tree1 = ParseNewick("(((A,B),C),(D,E));")
    tree2 = ParseNewick("((A,C),(B,D,E));")
    tree3 = ParseNewick("(((B,C),A),D,E);")
    trees = [tree1, tree2, tree3]
    
    result = newick(ParseNewick("((A,B,C),D,E);"))
    @test newick(majority_consensus_tree(trees)) == result

    tree4 = ParseNewick("(((B,C),A),D,F);")
    push!(trees, tree4)
    @test_throws ArgumentError majority_consensus_tree(trees)

    tree1 = ParseNewick("(((A,B)AB,C)ABC,(D,E)DE)no_name;")
    tree2 = ParseNewick("((A,C)AC,(B,D,E)BDE)no_name;")
    tree3 = ParseNewick("(((B,C)BC,A)BCA,D,E)no_name;")
    trees = [tree1, tree2, tree3]
    
    result = newick(ParseNewick("((A,B,C),D,E);"))
    @test newick(majority_consensus_tree(trees)) == result
end

@testset "loose_consensus_tree" begin
    tree1 = ParseNewick("(((A,B),C),(D,E));")
    tree2 = ParseNewick("((A,C),(B,D,E));")
    tree3 = ParseNewick("(((B,C),A),D,E);")
    trees = [tree1, tree2, tree3]
    
    result = newick(ParseNewick("(A,B,C,(D,E));"))
    @test newick(loose_consensus_tree(trees)) == result

    tree4 = ParseNewick("(((B,C),A),D,F);")
    push!(trees, tree4)
    @test_throws ArgumentError loose_consensus_tree(trees)
    tree1 = ParseNewick("(((A,B)AB,C)ABC,(D,E)DE)no_name;")
    tree2 = ParseNewick("((A,C)AC,(B,D,E)BDE)no_name;")
    tree3 = ParseNewick("(((B,C)BC,A)BCA,D,E)no_name;")
    trees = [tree1, tree2, tree3]
    
    result = newick(ParseNewick("(A,B,C,(D,E)no_name);"))
    @test newick(loose_consensus_tree(trees)) == result
end

@testset "greedy_consensus_tree" begin
    tree1 = ParseNewick("(((A,B),C),(D,E));")
    tree2 = ParseNewick("((A,C),(B,D,E));")
    tree3 = ParseNewick("(((B,C),A),D,E);")
    trees = [tree1, tree2, tree3]
    
    result = newick(ParseNewick("(((A,B),C),(D,E));"))
    # result = newick(ParseNewick("(((A,C),B),(D,E));"))
    # result = newick(ParseNewick("(((B,C),A),(D,E));"))
    @test newick(greedy_consensus_tree(trees)) == result
end

@testset "count_cluster_occurences" begin
    bit_vec = BitArray.([[1,0,0,0,0], [0,1,0,0,0], [0,0,1,0,0], [0,0,0,1,0],
              [0,0,0,0,1], [1,1,0,0,0], [1,1,1,0,0], [0,0,0,1,1], [1,0,0,0,0],
              [0,1,0,0,0], [0,0,1,0,0], [0,0,0,1,0], [0,0,0,0,1], [1,0,1,0,0],
              [0,1,0,1,1], [0,1,1,0,0], [1,1,1,0,0], [0,0,0,1,1], [1,0,0,0,0],
              [0,1,0,0,0], [0,0,1,0,0], [0,0,0,1,0], [0,0,0,0,1]])
    sort!(bit_vec, alg=QuickSort)
    queue = MCPhyloTree.PriorityQueue{BitVector, Int64}()
    queue[BitArray([0,0,0,0,1])] = -3
    queue[BitArray([0,0,0,1,0])] = -3
    queue[BitArray([0,1,0,0,0])] = -3
    queue[BitArray([0,0,1,0,0])] = -3
    queue[BitArray([1,0,0,0,0])] = -3
    queue[BitArray([0,0,0,1,1])] = -2
    queue[BitArray([1,1,0,0,0])] = -1
    queue[BitArray([1,0,1,0,0])] = -1
    queue[BitArray([0,1,0,1,1])] = -1
    queue[BitArray([0,1,1,0,0])] = -1
    queue[BitArray([1,1,1,0,0])] = -2

    @test MCPhyloTree.count_cluster_occurences(bit_vec) == queue
end
