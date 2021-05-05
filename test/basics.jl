
@testset "insert_node!" begin
    tree = MCPhyloTree.parsing_newick_string("(A,B,(C,D,E)F)G;")
    insert_tree = MCPhyloTree.parsing_newick_string("(A,B,((C,D,E)no_name)F)G;")
    insert_tree_newick = newick(insert_tree)
    insert_tree2 = MCPhyloTree.parsing_newick_string("(A,B,(((C,D)no_name,E)no_name)F)G;")
    insert_tree_newick2 = newick(insert_tree2)
    MCPhyloTree.number_nodes!(tree)

    children = [find_by_name(tree, "C"), find_by_name(tree, "D"),
                find_by_name(tree, "E")]
    mother = find_by_name(tree, "F")
    inserted_node = insert_node!(mother, children)
    @test newick(tree) == insert_tree_newick

    pop!(children)
    insert_node!(inserted_node, children)
    @test newick(tree) == insert_tree_newick2

    new_node = Node("G")
    append!(children, new_node)
    @test_throws AssertionError insert_node!(mother, children)
 end

@testset "delete_node!" begin
    tree = MCPhyloTree.parsing_newick_string("(A,B,(C,D)E)F;")
    MCPhyloTree.number_nodes!(tree)
    remove_tree = MCPhyloTree.parsing_newick_string("(A,B,C,D)F;")
    remove_tree_newick = newick(remove_tree)
    remove_tree2 = MCPhyloTree.parsing_newick_string("(B,C,D)F;")
    remove_tree2_newick = newick(remove_tree2)

    delete_node!(find_by_name(tree, "E"))
    @test newick(tree) == remove_tree_newick
    delete_node!(find_by_name(tree, "A"))
    @test newick(tree) == remove_tree2_newick
    @test_throws ArgumentError delete_node!(find_by_name(tree, "F"))
end

@testset "find_lca" begin
    tree = MCPhyloTree.parsing_newick_string("(A,B,(C,(D,E)F)G)H;")
    MCPhyloTree.set_binary!(tree)

    @testset "find_lca by name" begin
        @test find_lca(tree, ["A","C"]) == tree
        @test find_lca(tree, ["D", "E"]) == find_by_name(tree, "F")
        @test find_lca(tree, ["D", "F"]) == find_by_name(tree, "F")
        @test find_lca(tree, ["C", "D"]) == find_by_name(tree, "G")
        @test find_lca(tree, ["A","B","C"]) == tree
        @test find_lca(tree, ["D", "E", "F", "G"]) == find_by_name(tree, "G")
        @test find_lca(tree, ["D", "E", "C"]) == find_by_name(tree, "G")
    end # find_lca by name

    @testset "find_lca by nodes" begin
        nodes = [find_by_name(tree, "A"), find_by_name(tree, "B"),
                 find_by_name(tree, "C")]
        @test find_lca(tree, nodes) == tree

        nodes = [find_by_name(tree, "D"), find_by_name(tree, "E"),
                 find_by_name(tree, "F"), find_by_name(tree, "G")]
        @test find_lca(tree, nodes) == find_by_name(tree, "G")
    end # find_lca by nodes
end

@testset "check_leafsets" begin
    tree1 = MCPhyloTree.parsing_newick_string("(((A,B),C),(D,E))")
    tree2 = MCPhyloTree.parsing_newick_string("((A,C),(B,D,E))")
    tree3 = MCPhyloTree.parsing_newick_string("(((G,C),A),D,F)")
    tree4 = MCPhyloTree.parsing_newick_string("(((B,C),A),D,E)")
    tree5 = MCPhyloTree.parsing_newick_string("(((B,C),A),D,F)")

    trees = [tree1, tree2, tree3, tree4, tree5]
    MCPhyloTree.number_nodes!.(trees)
    MCPhyloTree.set_binary!.(trees)

    @test_throws ArgumentError MCPhyloTree.check_leafsets(trees)
end

@testset "mother" begin
    tree1 = MCPhyloTree.parsing_newick_string("(((A,B)F,C)G,(D,E)H)R;")
    MCPhyloTree.number_nodes!(tree1)
    MCPhyloTree.set_binary!(tree1)
    @test_throws ArgumentError MCPhyloTree.get_mother(find_by_name(tree1, "R"))
    @test MCPhyloTree.get_mother(find_by_name(tree1, "A")).name == "F"
end

