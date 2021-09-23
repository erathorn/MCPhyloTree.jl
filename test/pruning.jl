tree = ParseNewick("(A,B,(C,(D,E)F)G)H;")

@testset "Pruning leaf" begin
    node_list = ["C"]
    pt = prune_tree(tree, node_list)
    @test newick(pt) == newick(ParseNewick("(A,B,((D,E)F)G)H;"))

    node_list = ["C", "C"]
    @test_throws ArgumentError prune_tree(tree, node_list)

    node_list = ["K", "L"]
    @test_throws ArgumentError prune_tree(tree, node_list)

    node_list = ["L", "C"]
    pt = prune_tree(tree, node_list)
    @test newick(pt) == newick(ParseNewick("(A,B,((D,E)F)G)H;"))
end

@testset "Pruning root" begin
    node_list = ["H"]
    @test_throws ArgumentError prune_tree(tree, node_list)
end

@testset "Pruning misc" begin
    node_list1 = ["B", "C"]
    node_list2 = ["F"]
    node_list3 = ["A", "C", "G"]

    prune_tree1 = prune_tree(tree, node_list1)
    prune_tree2 = prune_tree(tree, node_list2)
    prune_tree3 = prune_tree(tree, node_list3)

    @test newick(prune_tree1) == newick(ParseNewick("(A,((D,E)F)G)H;"))
    @test newick(prune_tree2) == newick(ParseNewick("(A,B,(C)G)H;"))
    @test newick(prune_tree3) == newick(ParseNewick("(B)H;"))
end

@testset "Pruning inplace" begin
    node_list1 = ["B", "C"]
    prune_tree!(tree, node_list1)
    @test newick(tree) ==  newick(ParseNewick("(A,((D,E)F)G)H;"))
end
