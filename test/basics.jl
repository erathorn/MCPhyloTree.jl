
@testset "add_child!" begin
    tree = MCPhyloTree.parsing_newick_string("(A,B,(C,D,E)F)G;")
    MCPhyloTree.number_nodes!(tree)
    to_add = MCPhyloTree.parsing_newick_string("no_name;")
    println(typeof(to_add))
    root = MCPhyloTree.find_by_name(tree,"G")
    MCPhyloTree.add_child!(root,to_add)
    @test root.children[end].name == "no_name"
    @test length(MCPhyloTree.post_order(tree))==8

    tree = MCPhyloTree.parsing_newick_string("(A,B,(C,D,E)F)G;")
    MCPhyloTree.number_nodes!(tree)
    to_add = MCPhyloTree.parsing_newick_string("no_name;")
    println(typeof(to_add))
    leaf = MCPhyloTree.find_by_name(tree,"D")
    MCPhyloTree.add_child!(leaf,to_add)
    @test leaf.children[end].name == "no_name"
    @test length(MCPhyloTree.post_order(tree))==8


    tree = MCPhyloTree.parsing_newick_string("(A,B,(C,D,E)F)G;")
    MCPhyloTree.number_nodes!(tree)
    to_add = MCPhyloTree.parsing_newick_string("no_name;")
    println(typeof(to_add))
    mid = MCPhyloTree.find_by_name(tree,"F")
    MCPhyloTree.add_child!(mid,to_add)

    @test mid.children[end].name == "no_name"
    @test length(MCPhyloTree.post_order(tree))==8
end

@testset "remove_child!" begin
    tree = MCPhyloTree.parsing_newick_string("((A,B)C,(D,E)F)G;")
    MCPhyloTree.number_nodes!(tree)
    mother = MCPhyloTree.find_by_name(tree,"C")
    MCPhyloTree.remove_child!(mother,true)
    @test length(mother.children) == 1
    @test mother.children[1].name == "B"
    @test length(MCPhyloTree.post_order(tree)) == 6

    tree = MCPhyloTree.parsing_newick_string("((A,B)C,(D,E)F)G;")
    MCPhyloTree.number_nodes!(tree)
    mother = MCPhyloTree.find_by_name(tree,"C")
    MCPhyloTree.remove_child!(mother,false)
    @test length(mother.children) == 1
    @test mother.children[1].name == "A"
    @test length(MCPhyloTree.post_order(tree)) == 6

    tree = MCPhyloTree.parsing_newick_string("((A,B)C,(D,E)F)G;")
    MCPhyloTree.number_nodes!(tree)
    MCPhyloTree.remove_child!(tree,true)
    @test tree.children[1].name == "F"
    @test length(tree.children) == 1
    @test length(MCPhyloTree.post_order(tree)) == 4

end

@testset "remove_child! by node" begin
    tree = MCPhyloTree.parsing_newick_string("((A,B)C,(D,E)F)G;")
    MCPhyloTree.number_nodes!(tree)
    mother = MCPhyloTree.find_by_name(tree,"C")
    to_remove = MCPhyloTree.find_by_name(tree,"A")
    MCPhyloTree.remove_child!(mother,to_remove)
    @test length(mother.children) == 1
    @test mother.children[1].name == "B"
    @test length(MCPhyloTree.post_order(tree)) == 6

    tree = MCPhyloTree.parsing_newick_string("((A,B)C,(D,E)F)G;")
    MCPhyloTree.number_nodes!(tree)
    mother = MCPhyloTree.find_by_name(tree,"C")
    to_remove = MCPhyloTree.find_by_name(tree,"B")
    MCPhyloTree.remove_child!(mother,to_remove)
    @test length(mother.children) == 1
    @test mother.children[1].name == "A"
    @test length(MCPhyloTree.post_order(tree)) == 6

    tree = MCPhyloTree.parsing_newick_string("((A,B)C,(D,E)F)G;")
    MCPhyloTree.number_nodes!(tree)
    to_remove = MCPhyloTree.find_by_name(tree,"C")
    MCPhyloTree.remove_child!(tree,to_remove)
    @test tree.children[1].name == "F"
    @test length(tree.children) == 1
    @test length(MCPhyloTree.post_order(tree)) == 4
end

@testset "tree_length" begin
    tree = MCPhyloTree.parsing_newick_string("((A:5,B:5)C:3.5,(D:5,E:5)F:5)G:5;")
    @test MCPhyloTree.tree_length(tree) == 33.5
end

@testset "tree_height" begin
    tree = MCPhyloTree.parsing_newick_string("((A:5,B:5)C:3.5,(D:5,E:5)F:5)G:5;")
    @test MCPhyloTree.tree_height(tree) == 7
end

@testset "node_height" begin
    tree = MCPhyloTree.parsing_newick_string("((A:5,B:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.number_nodes!(tree)
    println("test")
    MCPhyloTree.node_height(tree)
    @test tree.height == 14.0
end

@testset "node_depth" begin
    tree = MCPhyloTree.parsing_newick_string("((A:5,B:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.number_nodes!(tree)
    MCPhyloTree.set_binary!(tree)
    mid = MCPhyloTree.find_by_name(tree,"F")
    leaf = MCPhyloTree.find_by_name(tree, "A")
    @test MCPhyloTree.node_depth(tree) == 0
    @test MCPhyloTree.node_depth(mid) == 1
    @test MCPhyloTree.node_depth(leaf) == 2
end

@testset "get_path" begin
    tree = MCPhyloTree.parsing_newick_string("((A:5,B:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.number_nodes!(tree)
    mid = MCPhyloTree.find_by_name(tree,"F")
    leaf = MCPhyloTree.find_by_name(tree,"A")
    @test MCPhyloTree.get_path(tree,mid) == [6]
    @test MCPhyloTree.get_path(tree,leaf) == [1, 5]
end

@testset "path_length" begin
    tree = MCPhyloTree.parsing_newick_string("((A:5,B:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.number_nodes!(tree)
    mid = MCPhyloTree.find_by_name(tree,"F")
    leaf = MCPhyloTree.find_by_name(tree,"A")
    @test MCPhyloTree.path_length(tree,mid) == 5.0
    @test MCPhyloTree.path_length(tree,leaf) == 14.0
end

@testset "get_sister" begin
    tree = MCPhyloTree.parsing_newick_string("((A:5,B:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.number_nodes!(tree)
    MCPhyloTree.set_binary!(tree)
    mid = MCPhyloTree.find_by_name(tree,"F")
    leaf = MCPhyloTree.find_by_name(tree,"A")
    leafresult = MCPhyloTree.get_sister(leaf)
    @test leafresult.name = "B"
    midresult = MCPhyloTree.get_sister(mid)
    @test midresult.name = "C"
end

@testset "set_binary!" begin
    tree = MCPhyloTree.parsing_newick_string("((A:5,B:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.number_nodes!(tree)
    MCPhyloTree.set_binary!(tree)
    mid = MCPhyloTree.find_by_name(tree,"F")
    leaf = MCPhyloTree.find_by_name(tree,"A")
    @test tree.binary == "1"
    @test leaf.binary =="1,0,0"
    @test mid.binary == "1,1"
end

@testset "number_nodes!" begin
    tree = MCPhyloTree.parsing_newick_string("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.number_nodes!(tree)
    MCPhyloTree.set_binary!(tree)
    mid = MCPhyloTree.find_by_name(tree,"F")
    leaf = MCPhyloTree.find_by_name(tree,"A")
    @test tree.num == 7
    @test mid.num == 6
    @test leaf.num == 1
end

@testset "blv" begin
    tree = MCPhyloTree.parsing_newick_string("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.number_nodes!(tree)
    MCPhyloTree.set_binary!(tree)
    @test MCPhyloTree.get_branchlength_vector(tree) == [5.0, 5.0, 5.0, 5.0, 9.0, 5.0]
    MCPhyloTree.set_branchlength_vector!(tree,[1.0,1.0,1.0,1.0,1.0,1.0])
    @test MCPhyloTree.get_branchlength_vector(tree) == [1.0,1.0,1.0,1.0,1.0,1.0]
    @test MCPhyloTree.get_sum_seperate_length!(tree) == [2.0, 4.0, 0.0, 0.0]
end

@testset "internal_external" begin
    tree = MCPhyloTree.parsing_newick_string("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.number_nodes!(tree)
    MCPhyloTree.set_binary!(tree)
    map = MCPhyloTree.internal_external_map(tree)
    @test map == [0, 0, 0, 0, 1, 1]
    @test MCPhyloTree.internal_external(tree) == [0, 0, 0, 0, 1, 1]
end

@testset "check_binary" begin
    tree = MCPhyloTree.parsing_newick_string("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    @test MCPhyloTree.check_binary(tree)
    tree = MCPhyloTree.parsing_newick_string("((B:5)C:9,(D:5,E:5)F:5)G:5;")
    @test !MCPhyloTree.check_binary(tree)
    tree = MCPhyloTree.parsing_newick_string("((B:5,A:5,X:5)C:9,(D:5,E:5)F:5)G:5;")
    @test !MCPhyloTree.check_binary(tree)
end




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
