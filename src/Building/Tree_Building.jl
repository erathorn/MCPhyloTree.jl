
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

    initialize_tree!(root)

    return root
end # function create_tree_from_leaves



"""
    from_df(df::Array{Float64,2}, name_list::Vector{String})::FNode

This function takes an adjacency matrix and a vector of names
and turns it into a tree. No checks are performed.

Returns the root node of the tree.

* `df` : matrix with edge weights
* `name_list` : a list of names such that they match the column indices of the matrix
"""
function from_df(df::Array{Float64,2}, name_list::Vector{String})::FNode

    node_list::Vector{FNode} = [Node(String(i)) for i in name_list]

    for (col_index, col) in enumerate(eachcol(df))
        for (row_index, entry) in enumerate(col)
            if entry != 0
                # there is a branch connecting these nodes
                node_list[col_index].inc_length = entry
                add_child!(node_list[row_index], node_list[col_index])
            end # end if
        end # end for
    end # end for
    i::Int = 0
    # find the root
    for (ind, n) in enumerate(node_list)
        if n.root == true
            i = ind
            break
        end # end if
    end # for
    node::FNode = node_list[i]

    # do some bookeeping here and set the binary representation of the nodes
    initialize_tree!(node)
    return node
end # function from_df


struct division
    arr::Array{<:Real, 2}
    names::Vector{<:AbstractString}
    number::Vector{Int}
end

Base.length(divi::division) = length(divi.names)

function single_inds(x::Array{Tuple{R, R}})::Vector{R} where R<:Real
    y = R[]
    for i in 1:length(x)
        if x[i][1] ∉ y && x[i][2] ∉ y
            push!(y, x[i][1])
        end
    end
    y
end

function create_division(covmat, names, numbers, tol)::Tuple{division, division}
    true_inds = single_inds(getfield.(findall(covmat .<= tol), :I))
    false_inds = [i for i in 1:size(covmat, 1) if i ∉ true_inds]
    @assert length(true_inds) > 0
    @assert length(false_inds) > 0
    d1 = division(covmat[true_inds, true_inds], names[true_inds], numbers[true_inds])
    d2 = division(covmat[false_inds, false_inds], names[false_inds], numbers[false_inds])
    d1, d2
end

function div2node(d1::division)
    n1 = Node(d1.names[1], d1.arr[1])
    n1.num = d1.number[1]
    n1
end


"""
    function cov2tree(covmat::Array{<:Real, 2}, names::Vector{<:AbstractString}, numbers::Vector{Int64}; tol::Real=1e-7)::FNode

This function reconstructs a tree from a covariance matrix. It takes a covariance matrix, 
a vector of leaf names and a vector of node numbers as mandatory arguments. The order of the two vectors must
correspond to the order of rows and columns in the covariance matrix. Optionally, the `tol` paramter indicates 
the boundary below which all values are treated as zero.

Returns the root node of the tree corresponding to the supplied covariance matrix.

* `covmat` : covariance matrix
* `names` : a list of names such that they match the column/row indices of the matrix
* `numbers` : a list of Integers such that they match the column/row indices of the matrix
* `tol` : cut off value below which all values are treated as zero
"""
function cov2tree(covmat::Array{<:Real, 2}, names::Vector{<:AbstractString}, numbers::Vector{Int64}; tol::Real=1e-7)::FNode

    if all(covmat .> tol)
        covmat = covmat .- minimum(covmat)
    elseif any(covmat .< 0.0)
        covmat[covmat .< 0] .= 0.0
    end
    try
        cov2tree_int(covmat, names, numbers, tol=tol)
    catch
        println(covmat)
    end
        
        
end

function cov2tree_int(covmat::Array{R, 2}, names::Vector{<:AbstractString}, numbers; tol::Real=1e-7)::GeneralNode{R, Int64} where R<:Real
    
    d1, d2 = create_division(covmat, names, numbers, tol)
    if length(d1) == 1 && length(d2) == 1
        # two singlular nodes left
        n1 = div2node(d1)
        n2 = div2node(d2)
        
    elseif length(d1) == 1 && length(d2) > 1
        # d1 is a leaf node, d2 is a cluster of nodes
        n1 = div2node(d1)
        min = minimum(d2.arr)
        n2 = cov2tree_int(d2.arr .- min, d2.names, d2.number, tol=tol)
        n2.inc_length = min

    elseif length(d1) > 1 && length(d2) == 1
        # d1 is a cluster of node2, d2 is a leaf node
        min = minimum(d1.arr)
        n1 = cov2tree_int(d1.arr .- min , d1.names, d1.number, tol=tol)
        n1.inc_length = min
        n2 = div2node(d2)
    
    else
        # d1 and d2 are a cluster of nodes
        min = minimum(d1.arr)
        n1 = cov2tree_int(d1.arr .- min, d1.names, d1.number, tol=tol)
        n1.inc_length = min
        min = minimum(d2.arr)
        n2 = cov2tree_int(d2.arr .- min, d2.names, d2.number, tol=tol)
        n2.inc_length = min
    end
    n = Node("root", one(R))
    add_child!(n, n1)
    add_child!(n, n2)
    return n
end