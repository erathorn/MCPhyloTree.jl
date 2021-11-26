# TODO: check trees are on the same leave set

"""
    RF(tree1::T, tree2::T)::Int64 where T <:GeneralNode

Calculate the Robinson-Foulds distance between the two trees.
In its current form the function assumes the trees have identical leave sets.

Returns result of algorithm as integer.

* `tree1` : tree used to determine RF distance.

* `tree2` : tree used to determine RF distance.
"""
function RF(tree1::T, tree2::T)::Int64 where T <: GeneralNode
    set_binary!(tree1)
    set_binary!(tree2)
    number_nodes!(tree1)
    number_nodes!(tree2)
    r, _ = RF_int(tree1, tree2)
    return r
end

function RF_int(tree1::T, tree2::T)::Tuple{Int64, Int64} where T <: GeneralNode
    bt3 = get_bipartitions(tree1)
    bt4 = get_bipartitions(tree2)
    length(bt3)+length(bt4) - 2* length(intersect(bt3, bt4)), length(bt3)+length(bt4)
end

"""
    RF_weighted(tree1::T, tree2::T)::Float64 where T <:GeneralNode

Calculate the weighted Robinson-Foulds distance between the two trees.
The raw Robinson-Foulds distance is weighted by the maximum distance between the trees.
In its current form the function assumes the trees have identical leave sets.

Returns result of algorithm as integer.

* `tree1` : tree used to determine RF distance.

* `tree2` : tree used to determine RF distance.
"""
function RF_weighted(tree1::T, tree2::T)::Float64 where T <: GeneralNode
    set_binary!(tree1)
    set_binary!(tree2)
    number_nodes!(tree1)
    number_nodes!(tree2)
    rf, l1 = RF_int(tree1, tree2)
    rf / l1
end

"""
    get_bipartitions(tree::T)::Vector{Tuple} where T <:GeneralNode

Get a vector of all bipartions of `tree`.

Returns a vector containing Tuples of sets representing the bipartions.
"""
function get_bipartitions(tree::T)::Vector{Tuple} where T <:GeneralNode
    po_vect= post_order(tree)[1:end-1]
    bt = Vector{Tuple}(undef, length(po_vect))
    all_leaves = [n.name for n in get_leaves(tree)]
    for ind in eachindex(po_vect)
        elem = po_vect[ind]::T
        outset = String[]
        inset = sort([i.name for i in get_leaves(elem)])
        outset = setdiff(all_leaves, inset)
        @inbounds bt[ind] = (join(sort(inset),","), join(sort(outset),","))
    end # for
    bt
end

"""
    get_bipartitions_as_bitvectors(tree::T)::Vector{BitVector} where T<:GeneralNode

Get all bipartions of `tree`.

Returns a vector containing BitVectors representing the splits.
"""
function get_bipartitions_as_bitvectors(tree::T)::Vector{BitVector} where T<:GeneralNode
    po_vect= post_order(tree)[1:end-1]
    bt = Vector{BitVector}(undef, length(po_vect))
    l::Int64 = length(get_leaves(tree))
    for (ind, node) in enumerate(po_vect)
        bit_vector::BitVector = falses(l)
        for leaf in get_leaves(node)
            bit_vector[leaf.num] = 1
        end # for
        @inbounds bt[ind] = bit_vector
    end # for
    bt
end

"""
    get_split(node::T, l::Int64)::BitVector where T<:GeneralNode

Get the split that you get by splitting the tree at `node`. Needs the number of tree leaves
as 2nd argument.

Returns a BitVector that represents the split.
"""
function get_split(node::T, l::Int64)::BitVector where T<:GeneralNode
    split::BitVector = falses(l)
    for leaf in get_leaves(node)
        split[leaf.num] = 1
    end # for
    split 
end

"""
    BHV_bounds(tree1::T, tree2::T)::Tuple{Float64, Float64} where T <:GeneralNode

This function calculates the lower and upper bounds of the geodesic in the
Billera-Holmes-Vogtman space.

Returns tuple of floats.
"""
function BHV_bounds(tree1::T, tree2::T)::Tuple{Float64, Float64} where T <:GeneralNode
    res_upper_1 = 0.0
    res_upper_2 = 0.0
    res_upper_3 = 0.0
    po::Vector{T} = post_order(tree1)
    for node in po[1:end-1]
        nom = find_num(tree2, node.num)
        if get_mother(node).num == get_mother(nom).num
            res_upper_3 += (node.inc_length-nom.inc_length)^2
        else
            res_upper_2  += nom.inc_length^2
            res_upper_1 += node.inc_length^2
        end # if

    end # for

    res_low::Float64 = res_upper_1+res_upper_2+res_upper_3
    res_high::Float64 = (sqrt(res_upper_2)+sqrt(res_upper_1))^2+res_upper_3

    sqrt(res_low), sqrt(res_high)
end

"""
helper function produces vector of tuples summarizing tree; tuple[1] = node number,
    tuple[2] = node name, tuple[3] = branchlength vector[num],
    tuple[4] = reference to leaf nodes under a given node
    comparable to input on line 250 of polymain.java
"""
function tree_summary(tree::T) where T <:GeneralNode
    blv = get_branchlength_vector(tree)
    ret = Vector{Tuple}(undef,length(blv))

    for i=1:length(blv)
        cur_node = find_num(tree,i)
        ret[i]=(i,cur_node.name,blv[i],get_leaves(cur_node))
    end #for
    ret
end
"""
helper function produces vector of common edges as well as their associated length;
each element is a tuple, containing reference to a node from tree1, a node from tree2,
and the length associated.
comparable to line 306 in PolyMain.java
"""
function common_edges(tree1::T,tree2::T) where T<:GeneralNode
    ret = []
    po::Vector{T} = post_order(tree1)
    for node in po[1:end-1]
        nom = find_num(tree2, node.num)
        if get_mother(node).num == get_mother(nom).num
            push!(ret,((node,nom,abs(node.inc_length-nom.inc_length))))
        end # if
    end # for
    ret
end
"""
returns list of references to non-common edges
comparable to lines 298-303 of PolyMain.java
"""
function non_common_edges(tree1::T,tree2::T) where T<:GeneralNode
    edgenumbers = [x=x[1].num for x in common_edges(tree1,tree2)]
    po_list = MCPhyloTree.post_order(tree1)
    nocommonlist = [x for x in po_list if !(x.num in edgenumbers)]
    return nocommonlist
end

mutable struct Geodesic
    ratioSequence::Vector{Float64}
    eLeafAttribs::Vector{Float64}
    fLeafAttribs::Vector{Float64}
    leafContributionSquared::Float64
    commonEdges::Vector{Tuple{FNode,FNode,Float64}}
    Geodesic(rs::Vector{Float64}, eLengths::Vector{Float64}, fLengths::Vector{Float64}, lcs) = 
        new(rs, eLengths, fLengths, 0.0, Tuple{FNode,FNode,Float64}[])
end # struct


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
    c_e_lengths::Vector{Tuple{Float64, Float64}} = get_common_edge_lengths([tree1, tree2], 
                                                                           common_edges, 
                                                                           length(leaves1))
end # geodesic

function get_common_edge_lengths(trees::Vector{FNode}, 
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
    leaves1::Vector{FNode} = sort!(get_leaves(tree1), by=x->x.num)
    leaves2::Vector{FNode} = sort!(get_leaves(tree2), by=x->x.num)
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

    common_edge::Tuple{FNode, Float64} = common_edges[1]
    split::BitVector = get_split(common_edge[1], length(leaves1))
    

    # find the common node in each tree
    common_node1, rev1 = get_node_from_split(tree1, split, leaves1)
    common_node2, rev2 = get_node_from_split(tree2, split, leaves2)

    # tracks if we have to swap the tree pairs after using the reversed split for one of 
    # the common nodes 
    reverse::Bool = rev1 ⊻ rev2
    
    # split the trees at the common nodes
    # tree1 = split_tree!(common_node1, tree1)
    # tree2 = split_tree!(common_node2, tree2)



    mother1 = common_node1.mother
    mother2 = common_node2.mother
    remove_child!(mother1, common_node1)
    remove_child!(mother2, common_node2)
    common_node1.root = true
    common_node2.root = true
    println([leaf.name for leaf in get_leaves(mother1)])
    println([leaf.name for leaf in get_leaves(mother2)])
    
    if mother1.root 
        if mother1.nchild == 1
            tree1 = remove_child!(mother1, mother1.children[1])
            tree1.root = true
        end # if
    else
        tree1 = MCPhyloTree.reroot(tree1, mother1)
        for child in tree1.children
            if child.nchild == 1
                inc_l = child.inc_length
                new_child = remove_child!(child, child.children[1])
                remove_child!(tree1, child)
                add_child!(tree1, new_child)
                new_child.inc_length += inc_l
            end # if /else
        end # for
        for leaf in get_leaves(tree1)
            if leaf.mother.nchild == 1 && !leaf.mother.root
                grandma = leaf.mother.mother
                inc_l = leaf.mother.inc_length
                remove_child!(grandma, leaf.mother)
                remove_child!(leaf.mother, leaf)
                add_child!(grandma, leaf)
                leaf.inc_length += inc_l
            end # if
        end # for
    end # if / else
    if mother2.root 
        if mother2.nchild == 1
            tree2 = remove_child!(mother2, mother2.children[1])
            tree2.root = true
        end # if 
    else
        tree2 = MCPhyloTree.reroot(tree2, mother2)
        for child in tree2.children
            if child.nchild == 1
                inc_l = child.inc_length
                new_child = remove_child!(child, child.children[1])
                new_child.inc_length += inc_l
                remove_child!(tree2, child)
                add_child!(tree2, new_child)
            end # if /else
        end # for
        for leaf in get_leaves(tree2)
            if leaf.mother.nchild == 1 && !leaf.mother.root
                grandma = leaf.mother.mother
                inc_l = leaf.mother.inc_length
                remove_child!(grandma, leaf.mother)
                remove_child!(leaf.mother, leaf)
                add_child!(grandma, leaf)
                leaf.inc_length += inc_l
            end # if
        end # for
    end # if / else
    # tree1 = mother1.root ? tree1 : MCPhyloTree.reroot(tree1, mother1) 
    # tree2 = mother2.root ? tree2 : MCPhyloTree.reroot(tree2, mother2) 
      
    
    # TODO: after rebase, maybe use new initialize method
    # set the nodes where we split the tree as the roots of the subtree
    # common_node1.root = true
    # common_node2.root = true
    set_binary!(common_node1)
    set_binary!(common_node2)
    set_binary!(tree1)
    set_binary!(tree2)
    number_nodes!(common_node1)
    number_nodes!(common_node2)
    number_nodes!(tree1)
    number_nodes!(tree2)

    fIO =  open("example.txt","a")
    if reverse
        write(fIO, newick(tree1))
        write(fIO, "\n")
        write(fIO, newick(common_node2))
        write(fIO, "\n")
        write(fIO, newick(tree2))
        write(fIO, "\n")
        write(fIO, newick(common_node1))
        write(fIO, "\n\n")
    else
        write(fIO, newick(tree1))
        write(fIO, "\n")
        write(fIO, newick(tree2))
        write(fIO, "\n")
        write(fIO, newick(common_node1))
        write(fIO, "\n")
        write(fIO, newick(common_node2))
        write(fIO, "\n\n")
    end # if / else

    close(fIO)
    println(length(get_leaves(common_node1)))
    println(length(get_leaves(common_node2)))
    println(length(get_leaves(tree1)))
    println(length(get_leaves(tree2)))

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
    split_tree!(split_node::FNode, tree::FNode)::FNode

This function splits a input tree at a node, and makes sure that internal nodes with 0 or 1
children after the split, are either deleted or fused with another node.

Returns the (potentially changed) root node of the subtree that is not headed by the node
where the tree was split.

* `split_node : The node where the tree will be split

* `tree` : root node of the tree
"""
function split_tree!(split_node::FNode, tree::FNode)::FNode
    mother::FNode = split_node.mother
    remove_child!(mother, split_node)
    if mother.nchild == 0 && !mother.root
        remove_child!(mother.mother, mother)
    end # if
    if mother.nchild == 1
        child = remove_child!(mother, mother.children[1])
        if mother.root
            child.root = true
            tree = child
        else
        add_child!(mother.mother, child)
        remove_child!(mother.mother, mother)
        end # if/else
    end # if
    tree
end # split_tree!


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
        elseif isCompatibleWith(split, tree_splits)
            length_diff = node.inc_length
            push!(commonEdges, (node, length_diff))
        end # elseif
    end # for

    tree_splits = get_bipartitions_as_bitvectors(tree1)
    for node in post_order(tree2)
        (node.nchild == 0 || node.root) && continue
        split::BitVector = get_split(node, l)
        if isCompatibleWith(split, tree_splits) && !(split in tree_splits)
            length_diff::Float64 = node.inc_length 
            push!(commonEdges, (node, length_diff))
        end # if
    end # for
    return commonEdges
end # getCommonEdges


"""
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


"""
    isCompatibleWith(b1::BitVector, b2::BitVector)::Bool

--- INTERNAL ---
Check if two splits are compatible.

Returns true if `b1` & `b2` - the bit vectors representing the two splits - are compatible.
"""
function isCompatibleWith(b1::BitVector, b2::BitVector)::Bool
    empty_intersect1::Bool = all(.!(b1 .& b2))
    empty_intersect2::Bool = all(.!(b1 .& .!b2))
    empty_intersect3::Bool = all(.!(.!b1 .& b2))
    empty_intersect4::Bool = all(b1 .| b2)
    return empty_intersect1 || empty_intersect2 || empty_intersect3 || empty_intersect4   
end # isCompatibleWith
