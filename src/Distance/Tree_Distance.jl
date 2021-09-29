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
    comparable to input on line 250 of polymain.java in gtp_170317\source code\polyAlg\PolyMain.java
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
comparable to line 306 of gtp_170317\source code\polyAlg\PolyMain.java
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
    leafVector1::Vector{String} = sort!([leaf.name for leaf in leaves1])
    leafVector2::Vector{String} = sort!([leaf.name for leaf in leaves2])
    leafVector1 != leafVector2 && 
        throw(ArgumentError("The two input trees do not have the same sets of leaves")) 
    for (leaf1, leaf2) in zip(leaves1, leaves2)
        abs(leaf1.inc_length - leaf2.inc_length) ^ 2
    end # for

    
    
    geo = Geodesic(RatioSequence(), t1LeafEdgeAttribs, t2LeafEdgeAttribs, 
                   leafContributionSquared, commonEdges)
end # geodesic