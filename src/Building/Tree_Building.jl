
function tree_from_leaves(leaf_nodes::Vector{String}, final_length::Int64)::Tuple{Vector{FNode}, Int}
    my_node_list::Array{FNode,1} = []

    # first create a list of leaf nodes
    for node_name in leaf_nodes
        nn = Node(node_name)

        push!(my_node_list,nn)
    end # for

    # Internal nodes are created using integers as names.
    temp_name::Int = length(my_node_list)+1

    # shuffle the node list to get a random tree
    Random.shuffle!(my_node_list)

    while length(my_node_list) > final_length
        # get two nodes
        # create a new mother node to which the two first nodes are added as children
        # add the new mother node to the list and reshuffle
        first_child::FNode = pop!(my_node_list)
        first_child.inc_length = rand()#*0.1
        second_child::FNode = pop!(my_node_list)
        second_child.inc_length = rand()
        curr_node::FNode = Node(string(temp_name))

        add_child!(curr_node, first_child)
        add_child!(curr_node, second_child)
        push!(my_node_list, curr_node)
        temp_name += 1
        Random.shuffle!(my_node_list)
    end # while

    return my_node_list, temp_name
end

"""
    function create_tree_from_leaves(leaf_nodes::Vector{String}, rooted::Bool=false)::FNode    
    
Build a random tree from a list of leaf names. The tree is unrooted by default.

Returns the root node of the new tree.

* `leaf_nodes` : A list of strings which are used as the names of the leaves.

* `rooted` : Boolean indicating if the tree should be rooted
"""
function create_tree_from_leaves(leaf_nodes::Vector{String}, rooted::Bool=false)::FNode
    
    my_node_list, temp_name = rooted ? tree_from_leaves(leaf_nodes, 2) : tree_from_leaves(leaf_nodes, 3)
    
    root::FNode = Node(string(temp_name))
    lchild = pop!(my_node_list)
    lchild.inc_length = rand()
    mchild = pop!(my_node_list)
    mchild.inc_length = rand()
    rchild = pop!(my_node_list)
    rchild.inc_length = rand()
    add_child!(root, lchild)
    add_child!(root, rchild)
    add_child!(root, mchild)

    set_binary!(root)
    number_nodes!(root)

    return root
end # function create_tree_from_leaves
