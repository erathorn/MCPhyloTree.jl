
@testset "Node Building" begin
    n = Node()
    @test isa(n, AbstractNode)    
    @test isa(n, GeneralNode{Float64, Int64})

    n = Node("Name")
    @test isa(n, AbstractNode)    
    @test isa(n, GeneralNode{Float64, Int64})
    @test n.name == "Name"

    n = Node("Name", 0.75)
    @test isa(n, AbstractNode)    
    @test isa(n, GeneralNode{Float64, Int64})
    @test n.name == "Name"
    @test n.inc_length == 0.75
end

@testset "Newick_Letters" begin
    tree = ParseNewick("(A,B)C;")
    @test length(tree) == 2
    @test tree.nchild == 2

    tree = ParseNewick("(A:1,B:1)C:1;")
    @test length(tree) == 2
    @test tree.nchild == 2
end

@testset "Newick_Numbers" begin
    tree = ParseNewick("(1,2)3;")
    @test length(tree) == 2
    @test tree.nchild == 2

    tree = ParseNewick("(1:1,2:1)3:1;")
    @test length(tree) == 2
    @test tree.nchild == 2
end

@testset "Newick_List" begin
    tree = ParseNewick("testtrees.nwk")
    @test length(tree) == 3

    @test_throws ArgumentError ParseNewick("I am no file")

    @test_throws ArgumentError ParseNewick("brokentree.nwk")
end

@testset "LargeNewick" begin
    Nleaves = 6000
    leave_list = ["leave_$i" for i in 1:Nleaves]

    gs_tree = MCPhyloTree.create_tree_from_leaves(leave_list, true)
    nstring = newick(gs_tree)
    tree = ParseNewick(nstring)
    
    @test gs_tree.height == tree.height
    @test length(collect(get_leaves(tree))) == Nleaves
    @test RF(gs_tree, tree) == 0
end



@testset "add_child!" begin
    tree = ParseNewick("(A,B,(C,D,E)F)G;")
    @test treesize(tree) == 7
    to_add = ParseNewick("no_name;")
    root = find_by_name(tree,"G")
    add_child!(root,to_add)
    @test root.children[end].name == "no_name"
    @test treesize(tree)==8

    tree = ParseNewick("(A,B,(C,D,E)F)G;")
    to_add = ParseNewick("no_name;")
    leaf = find_by_name(tree,"D")
    add_child!(leaf,to_add)
    @test leaf.children[end].name == "no_name"
    @test treesize(tree)==8


    tree = ParseNewick("(A,B,(C,D,E)F)G;")
    to_add = ParseNewick("no_name;")
    mid = find_by_name(tree,"F")
    add_child!(mid,to_add)

    @test mid.children[end].name == "no_name"
    @test treesize(tree)==8


    tree = ParseNewick("(A,B,(C,D,E)F)G;")
    to_add = ParseNewick("no_name;")
    mid = find_by_name(tree,"F")
    add_child!(mid,to_add, 1)

    @test mid.children[1].name == "no_name"
    @test treesize(tree)==8

    tree = ParseNewick("(A,B,(C,D,E)F)G;")
    to_add = ParseNewick("no_name;")
    mid = find_by_name(tree,"F")
    @test_throws AssertionError add_child!(mid,to_add, 5)

    tree = ParseNewick("(A,B,(C,D,E)F)G;")
    to_add = ParseNewick("no_name;")
    mid = find_by_name(tree,"F")
    
    add_child!(mid,to_add, 4)
    
    @test mid.children[4].name == "no_name"
    @test treesize(tree)==8

    tree = ParseNewick("(A,B,(C,D,E)F)G;")
    to_add = ParseNewick("no_name;")
    mid = find_by_name(tree,"F")
    add_child!(mid,to_add, 2)

    @test mid.children[2].name == "no_name"
    @test treesize(tree)==8



end

@testset "remove_child!" begin
    tree = ParseNewick("((A,B)C,(D,E)F)G;")
    mother = find_by_name(tree,"C")
    remove_child!(mother,true)
    @test length(mother.children) == 1
    @test mother.children[1].name == "B"
    @test treesize(tree) == 6

    tree = ParseNewick("((A,B)C,(D,E)F)G;")
    mother = find_by_name(tree,"C")
    remove_child!(mother,false)
    @test length(mother.children) == 1
    @test mother.children[1].name == "A"
    @test treesize(tree) == 6

    tree = ParseNewick("((A,B)C,(D,E)F)G;")
    remove_child!(tree,true)
    @test tree.children[1].name == "F"
    @test length(tree.children) == 1
    @test treesize(tree) == 4

end

@testset "remove_child! by node" begin
    tree = ParseNewick("((A,B)C,(D,E)F)G;")
    mother = find_by_name(tree,"C")
    to_remove = find_by_name(tree,"A")
    remove_child!(mother,to_remove)
    @test length(mother.children) == 1
    @test mother.children[1].name == "B"
    @test treesize(tree) == 6

    tree = ParseNewick("((A,B)C,(D,E)F)G;")
    mother = find_by_name(tree,"C")
    to_remove = find_by_name(tree,"B")
    remove_child!(mother,to_remove)
    @test length(mother.children) == 1
    @test mother.children[1].name == "A"
    @test treesize(tree) == 6

    tree = ParseNewick("((A,B)C,(D,E)F)G;")
    to_remove = find_by_name(tree,"C")
    remove_child!(tree,to_remove)
    @test tree.children[1].name == "F"
    @test length(tree.children) == 1
    @test treesize(tree) == 4
end

@testset "tree_length" begin
    tree = ParseNewick("((A:5,B:5)C:3.5,(D:5,E:5)F:5)G:5;")
    @test tree_length(tree) == 28.5
end

@testset "tree_height" begin
    tree = ParseNewick("((A:5,B:5)C:3.5,(D:5,E:5)F:5)G:5;")
    @test tree_height(tree) == 10.0
end

@testset "node_height" begin
    tree = ParseNewick("((A:5,B:5)C:9,(D:5,E:5)F:5)G:5;")
    node_height(tree)
    @test tree.height == 14.0
end

@testset "node_age" begin
    tree = ParseNewick("(((A:8,B:5)F:2,C:10)G:1,(D:12,E:4)H:3)R:1;")
    leaf = find_by_name(tree, "A")
    internal_node = find_by_name(tree, "F")
    furthest_leaf = find_by_name(tree, "D")
    @test node_age(furthest_leaf) == 0.0
    @test node_age(leaf) == 4.0
    @test node_age(internal_node) == 12.0
    @test node_age(tree) == 15.0
end 

@testset "node_depth" begin
    tree = ParseNewick("((A:5,B:5)C:9,(D:5,E:5)F:5)G:5;")
    mid = find_by_name(tree,"F")
    leaf = find_by_name(tree, "A")
    @test node_depth(tree) == 0
    @test node_depth(mid) == 1
    @test node_depth(leaf) == 2
end

@testset "get_path" begin
    tree = ParseNewick("((A:5,B:5)C:9,(D:5,E:5)F:5)G:5;")
    mid = find_by_name(tree,"F")
    leaf = find_by_name(tree,"A")
    @test get_path(tree,mid) == [6]
    @test get_path(tree,leaf) == [1, 5]
end

@testset "path_length" begin
    tree = ParseNewick("((A:5,B:5)C:9,(D:5,E:5)F:5)G:5;")
    mid = find_by_name(tree,"F")
    leaf = find_by_name(tree,"A")
    @test path_length(tree,mid) == 5.0
    @test path_length(tree,leaf) == 14.0
end

@testset "get_sister" begin
    tree = ParseNewick("((A:5,B:5)C:9,(D:5,E:5)F:5)G:5;")
    mid = find_by_name(tree,"F")
    leaf = find_by_name(tree,"A")
    leafresult = get_sister(leaf)
    @test leafresult.name == "B"
    midresult = get_sister(mid)
    @test midresult.name == "C"
end

@testset "set_binary!" begin
    tree = parsing_newick_string("((A:5,B:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.set_binary!(tree)
    mid = find_by_name(tree,"F")
    leaf = find_by_name(tree,"A")
    @test tree.binary == "1"
    @test leaf.binary =="1,0,0"
    @test mid.binary == "1,1"
end

@testset "number_nodes!" begin
    tree = parsing_newick_string("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.number_nodes!(tree)
    mid = find_by_name(tree,"F")
    leaf = find_by_name(tree,"A")
    @test tree.num == 7
    @test mid.num == 6
    @test leaf.num == 1
end

@testset "blv" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    @test get_branchlength_vector(tree) == [5.0, 5.0, 5.0, 5.0, 9.0, 5.0]
    set_branchlength_vector!(tree,[1.0,1.0,1.0,1.0,1.0,1.0])
    @test get_branchlength_vector(tree) == [1.0,1.0,1.0,1.0,1.0,1.0]
    @test MCPhyloTree.get_sum_seperate_length!(tree) == [2.0, 4.0, 0.0, 0.0]
end

@testset "internal_external" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    map = MCPhyloTree.internal_external_map(tree)
    @test map == [0, 0, 0, 0, 1, 1]
    @test MCPhyloTree.internal_external(tree) == [0, 0, 0, 0, 1, 1]
end

@testset "check_binary" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    @test check_binary(tree)
    tree = ParseNewick("((B:5)C:9,(D:5,E:5)F:5)G:5;")
    @test !check_binary(tree)
    tree = ParseNewick("((B:5,A:5,X:5)C:9,(D:5,E:5)F:5)G:5;")
    @test !check_binary(tree)
end

@testset "insert_node!" begin
    tree = ParseNewick("(A,B,(C,D,E)F)G;")
    insert_tree = ParseNewick("(A,B,((C,D,E)no_name)F)G;")
    insert_tree_newick = newick(insert_tree)
    insert_tree2 = ParseNewick("(A,B,(((C,D)no_name,E)no_name)F)G;")
    insert_tree_newick2 = newick(insert_tree2)

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
    tree = ParseNewick("(A,B,(C,D)E)F;")
    remove_tree = ParseNewick("(A,B,C,D)F;")
    remove_tree_newick = newick(remove_tree)
    remove_tree2 = ParseNewick("(B,C,D)F;")
    remove_tree2_newick = newick(remove_tree2)

    delete_node!(find_by_name(tree, "E"))
    @test newick(tree) == remove_tree_newick
    delete_node!(find_by_name(tree, "A"))
    @test newick(tree) == remove_tree2_newick
    @test_throws ArgumentError delete_node!(find_by_name(tree, "F"))
end

@testset "find_lca" begin
    tree = ParseNewick("(A,B,(C,(D,E)F)G)H;")

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
    tree1 = ParseNewick("(((A,B),C),(D,E));")
    tree2 = ParseNewick("((A,C),(B,D,E));")
    tree3 = ParseNewick("(((G,C),A),D,F);")
    tree4 = ParseNewick("(((B,C),A),D,E);")
    tree5 = ParseNewick("(((B,C),A),D,F);")

    trees = [tree1, tree2, tree3, tree4, tree5]

    @test_throws ArgumentError check_leafsets(trees)
end

@testset "mother" begin
    tree1 = ParseNewick("(((A,B)F,C)G,(D,E)H)R;")
    @test_throws ArgumentError get_mother(find_by_name(tree1, "R"))
    @test get_mother(find_by_name(tree1, "A")).name == "F"
end

@testset "force_ultrametric" begin
    tree1 = ParseNewick("(((A:5,B:5)F:5,C:5)G:5,(D:5,E:5)H:5)R:5;")
    MCPhyloTree.force_ultrametric!(tree1)
    for x in get_leaves(tree1)
        @test path_length(tree1,x) == 15.0
    end
end

@testset "ascii" begin
    tree = ParseNewick("(A,B)C;")
    singletree = ParseNewick("(((A,B)C)D,(E)F,G)H;")
    bigtree = ParseNewick("((A,B)C,(D,E)F)G;")
    lines,_ = MCPhyloTree.ascii(tree)
    biglines,_ = MCPhyloTree.ascii(bigtree)
    singlelines,_ = MCPhyloTree.ascii(singletree)
    @test lines == ["   /-A", "-C|", "   \\-B"]
    @test biglines == ["      /-A", "   /C|", "  |   \\-B", "-G|", "  |   /-D", "   \\F|", "      \\-E"]
    @test singlelines == ["         /-A", "   /D -C|", "  |      \\-B", "-H|", "  |-F --E", "  |", "   \\-G"]
end

@testset "childtype" begin
    test_node = Node()
    @test childtype(test_node) == GeneralNode{Float64, Int}
    @test childtype(typeof(test_node)) == GeneralNode{Float64, Int}

end