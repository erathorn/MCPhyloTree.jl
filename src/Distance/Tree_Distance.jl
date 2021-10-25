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
    all_leaves::Vector{String} = sort!([n.name for n in get_leaves(tree)])
    l::Int64 = length(all_leaves)
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
        end #if
    end #for
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
    t1LeafEdgeAttribs::Vector{Float64} = get_branchlength_vector(tree1) 
    t2LeafEdgeAttribs::Vector{Float64} = get_branchlength_vector(tree2)
    leaves1::Vector{FNode} = get_leaves(tree1)
    leaves2::Vector{FNode} = get_leaves(tree2)
    numEdges1::Int64 = length(t1LeafEdgeAttribs) - length(leaves1)
    numEdges2::Int64 = length(t2LeafEdgeAttribs) - length(leaves2)
    leafVector1::Vector{String} = sort!([leaf.name for leaf in leaves1])
    leafVector2::Vector{String} = sort!([leaf.name for leaf in leaves2])
    leafVector1 != leafVector2 && 
        throw(ArgumentError("The two input trees do not have the same sets of leaves")) 
    for (leaf1, leaf2) in zip(leaves1, leaves2)
        abs(leaf1.inc_length - leaf2.inc_length) ^ 2
    end # for

    
   """ 
    geo = Geodesic(RatioSequence(), t1LeafEdgeAttribs, t2LeafEdgeAttribs, 
                   leafContributionSquared, commonEdges)
    """
end # geodesic


function splitOnCommonEdge(tree1::FNode, tree2::FNode, numEdges1::Int64, numEdges2::Int64)
    numEdges1 == numEdges2 && return
    commonEdges::Vector{FNode} = getCommonEdges(tree1, tree2)
end # splitOnCommonEdge
    

function getCommonEdges(tree1::FNode, tree2::FNode)::Vector{Tuple{FNode, Float64}}
    commonEdges::Vector{Tuple{FNode, Float64}} = []
    # TODO maybe need to check leaves again; skipping for now 
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
            length_diff = node.inc_length - 
                                   find_lca(tree2, leaf_cluster).inc_length
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