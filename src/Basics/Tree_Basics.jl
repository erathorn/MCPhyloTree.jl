#=
my_tree:
- Julia version: 1.3.1
- Author: erathorn
- Date: 2019-05-07
=#

#TODO: Automate export of automatically genereated funtions

"""
    add_child!(mother_node::N, child::N, child_position::Union{Int64, Missing}=missing) where N <: GeneralNode

This function adds a child to the mother node.
The arity of the mother node is increased by `1` and the root
status of the child is set to `False`.

* `mother_node` : Node to add a child to.

* `child` : Node to add to mother_node.children.

* `child_position` : index at which to add the new child node; optional.
"""
function add_child!(mother_node::N, child::N, child_position::Union{Int64, Missing}=missing) where N <: GeneralNode
    if ismissing(child_position)
        push!(mother_node.children, child)
    else
        insert!(mother_node.children, child_position, child)
    end # if/else
    child.mother = mother_node
    mother_node.nchild += 1
    child.root = false
end # function add_child

"""
    remove_child!(mother_node::N, left::Bool)::N where N <: GeneralNode

This function removes a child from the list of nodes which are daughters of this
node. An input of "True" removes the left child,
while "False" removes the right child.

Returns the removed node.

* `mother_node` : Node from which to remove a child.

* `left` : boolean value determining which child of mother_node to remove.
"""
function remove_child!(mother_node::N, left::Bool)::N where N <: GeneralNode
    if left
        rv = popfirst!(mother_node.children)
        rv.mother = missing
    else
        rv = pop!(mother_node.children)
        rv.mother= missing
    end # end if
    mother_node.nchild -= 1
    return rv
end # function

"""
    remove_child!(mother_node::N, child::N)::N where N<:GeneralNode

This function removes a child from the list of nodes which are daughters of this
node.

The removed node is returned.

* `mother_node` : Node from which to remove a child.

* `child` : specific Node to remove. This node has to be a child of `mother_node`.
"""
function remove_child!(mother_node::N, child::N)::N where N<:GeneralNode
    ind = findfirst(x->x==child, mother_node.children)
    deleteat!(mother_node.children, ind)
    child.mother = missing

    mother_node.nchild -= 1
    return child
end # function

"""
    delete_node!(node::T)::Nothing where T<:GeneralNode

This functions deletes node from a tree and assigns all its children to its
mother node.

* `node` : Node to be deleted.
"""
function delete_node!(node::T)::Nothing where T<:GeneralNode
    if node.root == true
        throw(ArgumentError("Cannot remove root node"))
    end
    mother = get_mother(node)
    for child in node.children
        add_child!(mother, child, findfirst(x -> x == node, mother.children))
    end
    remove_child!(mother, node)
    return nothing
end

"""
    insert_node!(mother::T, children::Vector{T})::T where T<:GeneralNode

This function inserts a node into a tree after a mother node and gains
a subset of the mother's children as its children.

Returns the inserted node.

* `mother` : Node under which to add the newly-inserted node.

* `children` : Children of node referenced by "mother" to reassign as children of the newly-inserted node.
"""
function insert_node!(mother::T, children::Vector{T})::T where T<:GeneralNode
    @assert length(children) >= 1
    inserted_node = Node()
    for child in children
        @assert child in mother.children
    end # for
    index = findfirst(x -> x in children, mother.children)
    for child in children
        remove_child!(mother, child)
        add_child!(inserted_node, child)
    end # for
    add_child!(mother, inserted_node, index)
    return inserted_node
end




# legacy wrapper
function rescale_length(root::T) where T<:GeneralNode
    force_ultrametric(root)
end


"""
    force_ultrametric!(root::T) where T<:GeneralNode

Force an ultrametric version of the tree.
"""
function force_ultrametric!(root::T) where T<:GeneralNode
    po::Vector{T} = post_order(root)
    node2max_depth = zeros(UInt32, length(po))
    for node in po
        if node.nchild != 0
            mv = -1
            for child in node.children
                if node2max_depth[child.num] > mv
                    mv = node2max_depth[child.num]
                end # if
            end # for
            node2max_depth[node.num] = mv+1
        else
            node2max_depth[node.num] = 1
        end # if
    end # for

    node2dist = zeros(Float64, size(node2max_depth))
    tl = tree_height(root)
    nblv = zeros(Float64,length(po))
    for node in level_order(root)
        if node.root != true
            m = get_mother(node)
            nv = (tl - node2dist[m.num])/node2max_depth[node.num]

            nblv[node.num] = nv
            node2dist[node.num] = nv + node2dist[m.num]
        end # if
    end # for

    set_branchlength_vector!(root, nblv)
end # function force_ultrametric!


#################### Tree length & height ####################

"""
    tree_length(root::T)::Float64  where T<:GeneralNode

This function calculates the tree length.
"""
function tree_length(root::T)::Float64  where T<:GeneralNode
    return tree_length(root, 0.0)
end # function tree_length

"""
    tree_length(root::T, tl::Float64)::Float64 where T<:GeneralNode

This function does the internal tree length recursion
"""
function tree_length(root::T, tl::Float64)::Float64 where T<:GeneralNode

    #if length(root.children) != 0
    for child in root.children
        tl = tree_length(child, tl)
    end
    #end # if
    if root.root !== true
        tl += root.inc_length
    end

    tl
end # function tree_length


"""
    tree_height(root::T)::Float64  where T<:GeneralNode

This function calculates the tree height.
"""
function tree_height(root::T)::Float64  where T<:GeneralNode
    node_height(root)
    return root.height
end

"""
    node_height(root::T, mv::Float64)::Float64  where T<:GeneralNode

Calculate the height of a node.
"""
function node_height(root::T)  where T<:GeneralNode

    if root.nchild != 0
        for node in root.children
            node_height(node)
        end
        root.height = maximum([child.inc_length+child.height for child in root.children])
    else
        root.height = 0.0
    end
end # function node_height

function node_height_vec(root::T, vec::Vector{N})  where {T<:GeneralNode, N<:Real}

    if root.nchild != 0
        for node in root.children
            node_height_vec(node, vec)
        end
        root.height = maximum([child.inc_length+child.height for child in root.children])
    else
        root.height = 0.0
    end
    vec[root.num] = root.height
end # function node_height


function node_height_vec(root::T)::Vector{Float64} where T<:GeneralNode
    t = zeros(length(post_order(root)))
    node_height_vec(root, t)
    t
end # function node_height

"""
    function node_age(node::FNode)::Float64

Calculates the age of a node. If the tree is ultrametric then the node age is identical 
to the node height. It is calculated by subtracting the path length between the node & 
the root from the height of the root. This represents the age of the node, assuming the 
leaf farthest from the root has a node age of 0, and the root node is the 'oldest' node.
"""
function node_age(node::FNode)::Float64
    depth::Float64 = 0
    while !node.root
        depth += node.inc_length
        node = node.mother
    end
    return node.height - depth
end # function node_age

"""
    node_depth(node::T)::Int64 where T<:GeneralNode

Calculate the depth of a node.
"""
function node_depth(node::T)::Int64 where T<:GeneralNode
    return length(split(node.binary, ",")) - 1
end

function node_distance(tree::T, node1::T, node2::T)::Float64 where T<:GeneralNode
    lca = find_lca(tree, node1, node2)
    path_length(lca, node1)+path_length(lca,node2)
end

function get_path(ancestor::T, descendant::T)::Vector{Int64} where T<:GeneralNode
    path::Vector{Int64} = []
    while descendant.num != ancestor.num
        push!(path, descendant.num)
        descendant = get_mother(descendant)
    end
    path
end


"""
    path_length(ancestor::T, descendant::T)::Float64  where T<:GeneralNode

Note: The function assumes there is an ancestral relationship between the two nodes.

This function calculates the length of the path separating the ancestor from the
offspring node. The function follows the path specified through the binary
description of the node.
"""
function path_length(ancestor::T, descendant::T)::Float64  where T<:GeneralNode
    l::Float64 = 0.0

    while descendant != ancestor
        l += descendant.inc_length
        descendant = get_mother(descendant)
    end # while
    return l
end # function path_length


"""
    get_sister(node::T)::T  where T<:GeneralNode

This function gets the sister of `node`. It does so by looking for the respective
binary representation of the sister.
"""
@inline function get_sister(node::T)::T  where T<:GeneralNode
    m::T = get_mother(node)
    m.children[findfirst(y-> y!=node, m.children)]
end # function


"""
    get_mother(node::T)::T  where T<:GeneralNode

This function gets the mother of `node`. It does so by looking for the respective
binary representation of the mother node.
"""
@inline function get_mother(node::T)::T  where T <: GeneralNode
    ismissing(node.mother) && throw(ArgumentError("Node has no mother"))
    return node.mother
end # function

"""
    set_binary!(root::T)  where T <: GeneralNode

Assign a binary representation to each node, which specifies the path from the
root to this node via the binary representation of the node.
A left turn is a 1 in binary and a right turn a 0.
"""
function set_binary!(root::T)  where T <: GeneralNode
    if root.root
        root.binary = "1"
    end # if
    if root.nchild != 0
        for (ind, node) in enumerate(root.children)
            ind -= 1
            node.binary = string(root.binary,",", ind)
            set_binary!(node)
        end
    end # if
end # function set_binary

"""
    number_nodes!(root::T)::Nothing  where T<:GeneralNode

This function assigns a unique, sequential number to each node. Leaves are numbered first
in alphabetical order.
"""
function number_nodes!(root::T)::Nothing  where T<:GeneralNode
    tips = [n.name for n in get_leaves(root)]
    sort!(tips)
    running = length(tips)
    for node in post_order(root)
        if node.nchild == 0
            node.num = findfirst(x-> x==node.name, tips)
        else
            node.num = running + 1
            running += 1
        end
    end # for
end # function number_nodes

"""
    initialize_tree!(root::FNode; height::Bool=true)

This function initializes a tree, i.e. numbers its nodes and sets the binary + height 
fields.
"""
function initialize_tree!(root::FNode; height::Bool=true)
    set_binary!(root)
    number_nodes!(root)
    height && tree_height(root)
end # initialize_tree

"""
    update_tree!(root::FNode)

This function can be used to recompute the tree's binary and height values. This might be 
necessary after adding/moving/removing nodes.
"""
function update_tree!(root::FNode; height::Bool=true)
    set_binary!(root)
    height && tree_height(root)
end # update_tree

"""
    random_node(root::T)::T  where T<:GeneralNode

This function returns a random node from the tree.
"""
function random_node(root::T)::T  where T<:GeneralNode
    post_order_trav = post_order(root)
    return rand(post_order_trav)
end # function random_node



#################### Vector of branch lengths: get & set ####################

"""
    get_branchlength_vector(root::N)::Vector{T}  where {N <:GeneralNode, T<:Real}

Get the vector of branch lengths of the tree.
"""
function get_branchlength_vector(root::N)::Vector{Float64}  where {N <:GeneralNode}
    if length(root.blv) == 0
        root.blv = zeros(length(post_order(root))-1)
    end
    get_branchlength_vector(root, root.blv)
    return root.blv
end # function get_branchlength_vector


"""
    get_branchlength_vector(root::N, out_vec::Vector{T}) where {N<:GeneralNode, T<:Real}

Do post order traversal to retrieve a vector of branch lengths.
"""
function get_branchlength_vector(root::N, out_vec::Vector{T})::Nothing where {N<:GeneralNode, T<:Real}
    for child in root.children
        get_branchlength_vector(child, out_vec)
    end
    if !root.root
        out_vec[root.num] = root.inc_length
    end
    nothing
end



"""
    set_branchlength_vector!(root::N, blenvec::Array{T}) where {N<:GeneralNode, T<:Real}

This function sets the branch lengths of a tree to the values specified in blenvec.
"""
function set_branchlength_vector!(root::N, blenvec::Array{T}) where {N<:GeneralNode, T<:Real}
    any(0 .> blenvec) && throw("this should never happen")
    for child in root.children
        set_branchlength_vector!(child, blenvec)
    end

    @views if root.root !== true
        root.inc_length = blenvec[root.num]
    end
    nothing
end # function set_branchlength_vector!


"""
    get_sum_seperate_length!(root::T)::Vector{Float64}  where T<:GeneralNode

This function gets the sum of the branch lengths of the internal branches and the
branches leading to the leave nodes.
"""
function get_sum_seperate_length!(root::T)::Vector{Float64}  where T<:GeneralNode
    return get_sum_seperate_length!(post_order(root))
end # function get_sum_seperate_length!


"""
    get_sum_seperate_length!(post_order::Vector{T})::Vector{Float64}  where T<:GeneralNode

This function gets the sum of the branch lengths of the internal branches and the
branches leading to the leave nodes.
"""
function get_sum_seperate_length!(post_order::Vector{T})::Vector{Float64}  where T<:GeneralNode
    res_int::Float64 = 0.0
    res_leave::Float64 = 0.0
    res_int_log::Float64 = 0.0
    res_leave_log::Float64 = 0.0
    @simd for node in post_order
        if node.nchild !== 0
            # internal branches
            if !node.root
                res_int += node.inc_length
                res_int_log += log(node.inc_length)
            end
        else
            # branches leading to leaves
            res_leave += node.inc_length
            res_leave_log += log(node.inc_length)
        end # if
    end # for
    return [res_int, res_leave, res_int_log, res_leave_log]
end # function get_sum_seperate_length!

function internal_external_map(root::T)::Vector{Int64}  where T<:GeneralNode
    internal_external_map(post_order(root))
end

function internal_external_map(post_order::Vector{T})::Vector{Int64}  where T<:GeneralNode
    my_map::Vector{Int64} = zeros(Int64, length(post_order)-1)
    for node in post_order
        if !node.root
            if node.nchild != 0
                my_map[node.num] = 1
            end
        end
    end
    return my_map
end

function internal_external(root::T)::Vector{Int64}  where T<:GeneralNode
    v = root.IntExtMap
    length(v) != 0 && return v
    v = internal_external_map(root)
    root.IntExtMap = v

    return v
end


function find_lca(tree::T, node_l::Array{String, 1})::T  where T<:GeneralNode
    find_lca(tree, [find_by_name(tree, i) for i in node_l])
end

function find_lca(tree::T, node_l::Array{T})::T  where T<:GeneralNode
    @assert length(node_l) > 0
    if length(node_l) === 1
        return node_l[1]
    else
        n1 = popfirst!(node_l)
        n2 = popfirst!(node_l)
        lca = find_lca(tree, n1, n2)
        while length(node_l) !== 0
            n1 = popfirst!(node_l)
            lca = find_lca(tree, lca, n1)
        end
        return lca
    end
end

function find_lca(tree::T, node1::T, node2::T)::T  where T<:GeneralNode
    nb = lcp(node1.binary, node2.binary)
    thing = nb[end]
    if nb[end] == ','
        nb = nb[1:end-1]
    end
    if nb == node1.binary
        return node1
    end
    if nb == node2.binary
        return node2
    end
    find_binary(tree, nb)
end

"""
    check_binary(root::GeneralNode)::Bool

checks to see if given tree is binary; returns true if properly formatted and false otherwise
"""
function check_binary(root::GeneralNode)::Bool
    if root.root
        if root.nchild != 2 && root.nchild != 3
            return false
        end #if
    else
        if root.nchild != 0 && root.nchild != 2
            return false
        end #if
    end #else
    res::Bool = true
    for child in root.children
        res &= check_binary(child)
    end #for
    return res
end #function


"""
    check_leafsets(trees::Vector{T})::Nothing where T<:GeneralNode

Checks if an array of trees shares the same leafset (based on the leaf names)
"""
function check_leafsets(trees::Vector{T})::Nothing where T<:GeneralNode
    leaveset = Set([n.name for n in get_leaves(trees[1])])
    count = 0
    for (index, tree) in enumerate(trees[2:end])
        leaveset2 = Set([n.name for n in get_leaves(tree)])
        if leaveset != leaveset2
            println("Tree #$index has a different set of leaves than the first tree")
            count += 1
        end # if
    end # for
    if count != 0
        throw(ArgumentError("$count different trees have a different set of leaves than the first tree"))
    end # if
end # function

