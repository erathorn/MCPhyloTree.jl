# TODO: check trees are on the same leave set

"""
    RF(tree1::T, tree2::T)::Int64 where T <:AbstractNode

Calculate the Robinson-Foulds distance between the two trees.
In its current form the function assumes the trees have identical leave sets.

Returns result of algorithm as integer.

* `tree1` : tree used to determine RF distance.

* `tree2` : tree used to determine RF distance.
"""
function RF(tree1::T, tree2::T)::Int64 where T <: GeneralNode
    initialize_tree!(tree1, height=false)
    initialize_tree!(tree2, height=false)
    r, _ = RF_int(tree1, tree2)
    return r
end


function RF_int(tree1::T, tree2::T)::Tuple{Int64, Int64} where T <: GeneralNode
    bt3 = get_bipartitions(tree1)
    bt4 = get_bipartitions(tree2)
    length(bt3)+length(bt4) - 2* length(intersect(bt3, bt4)), length(bt3)+length(bt4)
end


"""
    RF_weighted(tree1::T, tree2::T)::Float64 where T <:AbstractNode

Calculate the weighted Robinson-Foulds distance between the two trees.
The raw Robinson-Foulds distance is weighted by the maximum distance between the trees.
In its current form the function assumes the trees have identical leave sets.

Returns result of algorithm as integer.

* `tree1` : tree used to determine RF distance.

* `tree2` : tree used to determine RF distance.
"""
function RF_weighted(tree1::T, tree2::T)::Float64 where T <: GeneralNode
    initialize_tree!(tree1, height=false)
    initialize_tree!(tree2, height=false)
    rf, l1 = RF_int(tree1, tree2)
    rf / l1
end


"""
    get_bipartitions(tree::T)::Vector{Tuple} where T <:AbstractNode

Get a vector of all bipartions of `tree`.

Returns a vector containing tuples of sets representing the bipartitions.
"""
function get_bipartitions(tree::T)::Vector{Tuple} where T <:AbstractNode
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

Get all bipartitions of `tree`.

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
    get_split(node::T, leaves_dict::Dict{String, Int64})::BitVector where T<:GeneralNode
    
This is a separate get_split function intended to ensure that that the bit vectors of 
subtrees can be compared.
    
Get the split that you get by splitting the tree at `node`. Needs a leaf dictionary
containing the name of the leaves and their number in a set tree, to ensure the bit vector 
can be compared to the bit vector of other subtrees.

Returns a BitVector that represents the split.
"""
function get_split(node::T, leaves2num::Dict{String, Int64})::BitVector where T<:GeneralNode         
    split::BitVector = falses(length(leaves2num))
    for leaf in get_leaves(node)
        split[leaves2num[leaf.name]] = 1
    end # for
    split 
end


"""
    BHV_bounds(tree1::T, tree2::T)::Tuple{Float64, Float64} where T <:AbstractNode

This function calculates the lower and upper bounds of the geodesic in the
Billera-Holmes-Vogtman space.

Returns tuple of floats.
"""
function BHV_bounds(tree1::T, tree2::T)::Tuple{Float64, Float64} where T <:AbstractNode
    bp_t1 = get_bipartitions(tree1)
    bp_t2 = get_bipartitions(tree2)
    
    T1minusT2 = 0.0
    T2minusT1 = 0.0
    T1andT2 = 0.0
    
    for bp in bp_t1
        found1 = find_lca(tree1, String.(split(bp[1], ",")))
        if found1.root
            found1 = find_lca(tree1, String.(split(bp[2], ",")))
        end
        
        ind = findfirst(isequal(bp), bp_t2)
        if !isnothing(ind)
            found2 = find_lca(tree2, String.(split(bp[1], ",")))
            if found2.root
                found2 = find_lca(tree2, String.(split(bp[2], ",")))
            end
            T1andT2 += (found1.inc_length - found2.inc_length)^2
            deleteat!(bp_t2, ind)
        else
            T1minusT2 += found1.inc_length^2
        end
    end
    
    # only nodes which are present in the second tree are left in bp_t2
    for bp in bp_t2
        found2 = find_lca(tree2, String.(split(bp[1], ",")))
        if found2.root
            found2 = find_lca(tree2, String.(split(bp[2], ",")))
        end
        T2minusT1 += found2.inc_length^2
    end
    
    
    res_low::Float64 = T1minusT2+T2minusT1+T1andT2
    res_high::Float64 = (sqrt(T1minusT2)+sqrt(T2minusT1))^2+T1andT2

    sqrt(res_low), sqrt(res_high)
end