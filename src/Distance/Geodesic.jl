const Ratio = Tuple{Vector{FNode}, Vector{FNode}}
const RatioSequence = Vector{Ratio}

mutable struct Geodesic
    ratioSequence::RatioSequence
    eLeafAttribs::Vector{Float64}
    fLeafAttribs::Vector{Float64}
    leafContributionSquared::Float64
    commonEdges::Vector{Tuple{FNode,FNode,Float64}}
    Geodesic(rs::RatioSequence, eLengths::Vector{Float64}, fLengths::Vector{Float64}) = 
        new(rs, eLengths, fLengths, 0.0, Tuple{FNode,FNode,Float64}[])

    Geodesic(rs::RatioSequence) = 
        new(rs, [], [], 0.0, Tuple{FNode,FNode,Float64}[])
end # Geodesic


function geodesic(tree1::FNode, tree2::FNode)
    leafContributionSquared::Float64 = 0.0
    leaves1::Vector{FNode} = get_leaves(tree1)
    leaves2::Vector{FNode} = get_leaves(tree2)
    leafVector1::Vector{String} = sort!([leaf.name for leaf in leaves1])
    leafVector2::Vector{String} = sort!([leaf.name for leaf in leaves2])
    leafVector1 != leafVector2 && 
        throw(ArgumentError("The two input trees do not have the same sets of leaves")) 
    for (leaf1, leaf2) in zip(leaves1, leaves2)
        leafContributionSquared += abs(leaf1.inc_length - leaf2.inc_length) ^ 2
    end # for
    
    """
    geo = Geodesic(RatioSequence(), t1LeafEdgeAttribs, t2LeafEdgeAttribs, 
                   leafContributionSquared, commonEdges)
    """
    non_common_edges::Vector{Tuple{FNode, FNode}} = splitOnCommonEdge(deepycopy(tree1), 
                                                                      deepcopy(tree2))

    common_edges::Vector{Tuple{FNode, Float64}} = getCommonEdges(tree1, tree2)
    c_e_lengths::Vector{Tuple{Float64, Float64}} = get_commonedge_lengths([tree1, tree2], 
                                                                           common_edges, 
                                                                           length(leaves1))
end # geodesic


function get_commonedge_lengths(trees::Vector{FNode}, 
                                common_edges::Vector{Tuple{FNode, Float64}}, l::Int64
                                )::Tuple{Vector{Float64}, Vector{Float64}}

    common_edge_lengths::Tuple{Vector{Float64}, Vector{Float64}} = ([], [])
    post_orders = post_order.(trees)
    bps = get_bipartitions_as_bitvectors.(trees)
    for common_edge in common_edges 
        split = get_split(common_edge[1], l)
        for i in [1,2]
            ind = findfirst(x -> x == split, bps[i])
            inc_length::Float64 = isnothing(ind) ? 0.0 : post_orders[i][ind].inc_length
            println(i)
            push!(common_edge_lengths[i], inc_length)
        end # for
    end # for
    common_edge_lengths
end # get_common_edge_lengths


"""
    splitOnCommonEdge(tree1::FNode, tree2::FNode; non_common_edges=[])
        ::Vector{Tuple{FNode, FNode}}

This function finds the common edges of two trees, and splits them at the first common edge
it finds. Then it recursively splits the resulting subtrees aswell. 

Returns a Vector of all pairs of subtrees that share no common edges.

* `tree1` : root node of the first tree

* `tree2` : root node of the second tree

* `non_common_edges` : vector of non common edges. It is initialized empty and appended to 
                       at each recursion
"""
function splitOnCommonEdge(tree1::FNode, tree2::FNode; non_common_edges=[]
                          )::Vector{Tuple{FNode, FNode}}

    numNodes1::Int64 = length(post_order(tree1))
    numNodes2::Int64 = length(post_order(tree2))
    leaves1::Vector{FNode} = sort!(get_leaves(tree1), by=x->x.name)
    leaves2::Vector{FNode} = sort!(get_leaves(tree2), by=x->x.name)
    numEdges1::Int64 = numNodes1 - length(leaves1) - 1
    numEdges2::Int64 = numNodes2 - length(leaves2) - 1
    (numEdges1 <= 0 || numEdges2 <= 0) && return []
    
    common_edges::Vector{Tuple{FNode, Float64}} = getCommonEdges(tree1, tree2)
    
    # if there are no common edges, add the trees to the array of subtrees that share no 
    # common edges
    if isempty(common_edges)
        push!(non_common_edges, (tree1, tree2))
        return non_common_edges
    end # if
    
    # get the first common edge that was found
    common_edge::Tuple{FNode, Float64} = common_edges[1]
    # get a bit vector representing the split of the common edge 
    split::BitVector = get_split(common_edge[1], length(leaves1))
    
    # find the common node in each tree by using its split
    common_node1, rev1 = get_node_from_split(tree1, split, leaves1)
    common_node2, rev2 = get_node_from_split(tree2, split, leaves2)

    # track if we have to swap the tree pairs after using the reversed split for one of 
    # the common nodes 
    reverse::Bool = rev1 ⊻ rev2
    
    # split the trees at the common nodes
    split_tree!(common_node1)
    split_tree!(common_node2)

    # initialize the new subtrees
    initialize_tree!(common_node1)
    initialize_tree!(common_node2)
    initialize_tree!(tree1)
    initialize_tree!(tree2)
    
    # recursively split new subtrees as well
    if reverse
        splitOnCommonEdge(tree1, common_node2; non_common_edges)
        splitOnCommonEdge(tree2, common_node1; non_common_edges)
    else
        splitOnCommonEdge(common_node1, common_node2; non_common_edges)
        splitOnCommonEdge(tree1, tree2; non_common_edges)
    end # if/else
    return non_common_edges
end # splitOnCommonEdge
 

"""
    getCommonEdges(tree1::FNode, tree2::FNode)::Vector{Tuple{FNode, Float64}}

This function returns all the common edges of two trees with the same leafset. It also
calculates the length difference of each pair of common edges. 

* `tree1` : root node of the first tree

* `tree2` : root node of the second tree
"""
function getCommonEdges(tree1::FNode, tree2::FNode)::Vector{Tuple{FNode, Float64}}
    commonEdges::Vector{Tuple{FNode, Float64}} = []
    # TODO maybe need to check leaves again; skipping for now 
    # TODO maybe return a split instead of a node
    tree_splits::Vector{BitVector} = get_bipartitions_as_bitvectors(tree2)
    l = length(get_leaves(tree1))
    leaves2 = get_leaves(tree2)
    for node in post_order(tree1)
        (node.nchild == 0 || node.root) && continue
        # get the split of the current node represented as a BitVector
        split::BitVector = get_split(node, l)
        length_diff::Float64 = 0.0
        if split in tree_splits
            leaf_cluster = leaves2[split]
            length_diff = node.inc_length - find_lca(tree2, leaf_cluster).inc_length
            push!(commonEdges, (node, length_diff))
        end
    end # for
    return commonEdges
end # getCommonEdges


"""
    get_node_from_split(tree::FNode, split::BitVector, leaves::Vector{FNode}
                       )::Tuple{FNode, Bool}

--- INTERNAL ---
This function finds the node in a tree that corresponds to the input split or the reverse of
the input split (in case the found node for the original split would have been the root).

Returns a Vector of all pairs of subtrees that share no common edges.

* `tree` : root node of the tree.

* `split` : A bit vector representing the split

* `leaves` : vector of the leaves of the tree
"""
function get_node_from_split(tree::FNode, split::BitVector, leaves::Vector{FNode}=[]
                            )::Tuple{FNode, Bool}

    if isempty(leaves)
        leaves = get_leaves(tree)
    end # if
    common_node::FNode = find_lca(tree, leaves[split])
    reverse::Bool = false
    # if the found node is the root, instead look for the lca of the other side of the split
    if common_node.root
        reverse = true
        # get the reversed bitsplit of the common edge, in case the common edge is the root
        common_node = find_lca(tree, leaves[(!).(split)])
    end # if
    (common_node, reverse)
end # get_node_from_split


"""
    split_tree!(node::FNode)::FNode

This function splits a tree at the input node. The input node acts as the root of one of the
resulting subtrees. The tree headed by the original root, receives a new node where the
split node used to be. It represents all the leaves that were below the split node.

* `node`: the node where the tree will be split.
"""
function split_tree!(node::FNode)::Nothing
    mother = node.mother
    # remove the split node from the tree
    remove_child!(mother, node)
    node.root = true
    # add node representing the leaves below the split node
    add_child!(mother, Node(join(sort!([leaf.name for leaf in get_leaves(node)]), " ")))
    return
end # split_tree!


"""
    get_geodesic_nocommonedges(tree1::FNode, tree2::FNode)
"""
function get_geodesic_nocommonedges(tree1::FNode, tree2::FNode)
    numNodes1::Int64 = length(post_order(tree1))
    numNodes2::Int64 = length(post_order(tree2))
    numEdges1::Int64 = numNodes1 - length(leaves1) - 1
    numEdges2::Int64 = numNodes2 - length(leaves2) - 1
    rs::RatioSequence = []
    aVertices::Vector{Int64} = []
    bVertices::Vector{Int64} = []
    queue::Vector{Ratio} = Vector{Ratio}()
    ratio::Ratio = ([],[])
    cover::Array{Int64} =[[]]

    commonedges = getCommonEdges(tree1, tree2)
    # doublecheck to make sure the trees have no common edges
    length(commonedges) == 0 && throw(ArgumentError("Exiting: Can't compute geodesic between subtrees that have common edges."))
    
    if numEdges1 == 0 || numEdges2 == 0
        throw(ArgumentError("Exiting: Can't compute the geodesic for trees with no edges"))
    end # if

    if numEdges1 == 1 || numEdges2 == 1
        internal_nodes1 = filter!(x -> x.nchild != 0, post_order(tree1)[1:end 1])
        internal_nodes2 = filter!(x -> x.nchild != 0, post_order(tree1)[1:end-1])
		push(rs, (internal_nodes1, internal_nodes2))
		return Geodesic(rs)
	end # if

    leaves::Vector{FNode} = get_leaves(tree1)
    leaf_dict::Dict{String, Int64} = Dict(leaf.name => leaf.num for leaf in leaves)
    internal_nodes1 = filter!(x -> x.nchild != 0, post_order(tree1)[1:end-1])
    internal_nodes2 = filter!(x -> x.nchild != 0, post_order(tree1)[1:end-1])
    incidence_matrix::BitMatrix =  get_incidence_matrix(internal_nodes1, internal_nodes2, 
                                                        leaf_dict)

    
end # get_geodesic_nocommonedges


"""
    get_incidence_matrix(edges1::Vector{FNode}, edges2::Vector{FNode}, 
        leaf_dict::Dict{String,Int64})::BitMatrix

--- INTERNAL ---
This function computes the incidence_matrix of two trees. A field of the matrix contains a 1
if the corresponding subtrees of our input trees are not compatible.

Returns the incidence_matrix.

* `edges1` : edges of the first tree

* `edges2` : edges of the second tree

* `leaf_dict` : Dictionary that links each leaf to a unique number
"""
function get_incidence_matrix(edges1::Vector{FNode}, edges2::Vector{FNode}, 
                              leaf_dict::Dict{String,Int64})::BitMatrix

    incidence_matrix::BitMatrix = trues((length(edges1), length(edges2)))   
    for i in 1:length(edges1)
        for j in 1:length(edges2)
            if crosses(get_split(edges1[i], leaf_dict), get_split(edges2[j], leaf_dict))
                incidence_matrix[i,j] = true
            else
                incidence_matrix[i,j] = false
            end # if/else
        end # for
    end # for
    return incidence_matrix
end # get_incidence_matrix


"""
    crosses(b1::BitVector, b2::BitVector)::Bool

--- INTERNAL ---
Check if two splits cross.

Returns true if `b1` & `b2` - the bit vectors representing the two splits - are neither
disjoint nor does either one contain the other.
"""
function crosses(b1::BitVector, b2::BitVector)::Bool
    disjoint::Bool = all(.!(b1 .& b2))
    contains1::Bool = all(b1 .>= b2)
    contains2::Bool = all(b2 .>= b1)
    return !(disjoint || contains1 || contains2)   
end # isCompatibleWith


#=
### might need these function at some point, probably not though

    isCompatibleWith(node_split::BitVector, bipartitions::Vector{BitVector})::Bool

--- INTERNAL ---
Check if a split is compatible with a vector of splits.

Returns true if `node_split` is compatible with all splits in `bipartitions`.
"""
function isCompatibleWith(node_split::BitVector, bipartitions::Vector{BitVector})::Bool
    for bipartition in bipartitions
        !isCompatibleWith(node_split, bipartition) && return false
    end # for
    true
end # isCompatibleWith
=#