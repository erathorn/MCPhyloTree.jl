mutable struct Ratio{T<:GeneralNode}
    e_length::Float64
	f_length::Float64
	e_edges::Vector{T}
	f_edges::Vector{T}

    Ratio() = new{GeneralNode}(0.0, 0.0, [], [])
    Ratio(e::Vector{T}, f::Vector{T}) where T<:GeneralNode = new{typeof(e[1])}(geo_avg(e), geo_avg(f), e, f)
end # Ratio


mutable struct RatioSequence
    ratios::Vector{Ratio}
    combine_code::Int64

    RatioSequence() = new([], 0)
end # RatioSequence


const EdgeLengths = Tuple{Float64, Float64}
const CommonEdge = Tuple{T, Float64} where T<:GeneralNode 


mutable struct Geodesic
    ratio_seq::RatioSequence
    e_leaf_attribs::Vector{Float64}
    f_leaf_attribs::Vector{Float64}
    leaf_contribution²::Float64
    common_edges::Vector{CommonEdge}
    common_edge_lengths::Vector{EdgeLengths}

    Geodesic(rs::RatioSequence, e_lengths::Vector{Float64}, f_lengths::Vector{Float64}) = 
        new(rs, e_lengths, f_lengths, 0.0, CommonEdge[], EdgeLengths[])

    Geodesic(rs::RatioSequence) = 
        new(rs, [], [], 0.0, CommonEdge[], EdgeLengths[])
end # Geodesic


function get_distance(geo::Geodesic)::Float64
    common_edge_dist²::Float64 = 0
    for i in 1:length(geo.common_edge_lengths)
        common_edge_dist² += (geo.common_edge_lengths[i][1] - 
                              geo.common_edge_lengths[i][2]) ^ 2
    end # for
    return sqrt(get_distance(get_non_des_rs_with_min_dist(geo.ratio_seq)) ^ 2 + common_edge_dist² + geo.leaf_contribution²)
end # get_distance


function get_distance(rs::RatioSequence)::Float64
    distance²::Float64 = 0
    r::Ratio = Ratio()
    for i in 1:length(rs.ratios)
        r = rs.ratios[i]
        distance² += (r.e_length + r.f_length) ^ 2
    end # for
    return sqrt(distance²)
end # get_distance


mutable struct Vertex
    label::Float64
    weight::Float64
    residual::Float64
    pred::Int64

    Vertex(weight::Float64) = new(0.0, weight*weight, 0.0, 0)
end # Vertex


mutable struct BipartiteGraph
    incidence_matrix::BitMatrix
	nums::Vector{Int64}
	a_vertices::Vector{Vertex}
    b_vertices::Vector{Vertex}
	debug::Bool
end # BipartiteGraph


function build_bipartite_graph(incidence_matrix::BitMatrix, a_weight::Vector{Float64}, 
                               b_weight::Vector{Float64})::BipartiteGraph

    a::Int64 = length(a_weight)
    b::Int64 = length(b_weight)
    n::Int64 = max(a, b)
    a_vertices::Vector{Vertex} = fill(Vertex(0.0), n)
    b_vertices::Vector{Vertex} = fill(Vertex(0.0), n)
    for i in 1:a
        a_vertices[i]= Vertex(a_weight[i])
    end # for 
    for i in 1:b
        b_vertices[i]= Vertex(b_weight[i])
    end # for
    return BipartiteGraph(incidence_matrix, [a, b, n, 0, 0], a_vertices, b_vertices, false)
end


function geodesic(tree1::T, tree2::T) where T<:GeneralNode
    trees::Vector{T} = [tree1, tree2]
    leaf_contribution²::Float64 = 0.0
    leaves::Vector{Vector{T}} = get_leaves.(trees)
    leaf_names::Vector{Vector{String}} = [[leaf.name for leaf in leaves[i]] for i in 1:2] 
    perms::Vector{Vector{Int64}} = sortperm.(leaf_names)
    leaf_vecs::Vector{Vector{String}} = [leaf_names[i][perms[i]] for i in 1:2]
    leaf_attribs::Vector{Vector{Float64}} = [[leaf.inc_length for leaf in leaves[i]][perms[i]] for i in 1:2]
    
    leaf_vecs[1] != leaf_vecs[2] && 
        throw(ArgumentError("The two input trees do not have the same sets of leaves")) 
    for (leaf1, leaf2) in zip(leaves[1], leaves[2])
        leaf_contribution² += abs(leaf1.inc_length - leaf2.inc_length) ^ 2
    end # for
    
    geo = Geodesic(RatioSequence(), leaf_attribs[1], leaf_attribs[2])
    geo.leaf_contribution² = leaf_contribution²

    non_common_edges::Vector{Tuple{T, T}} = split_on_common_edge(deepcopy.(trees)...)

    common_edges::Vector{CommonEdge} = get_common_edges(trees...)
    geo.common_edges = common_edges

    c_e_lengths::Vector{EdgeLengths} = get_common_edge_lengths(trees, common_edges,
                                                               length(leaves[1]))                                                     
    geo.common_edge_lengths = c_e_lengths

    for i in 1:length(non_common_edges)
        subtree_a = non_common_edges[i][1]
        subtree_b = non_common_edges[i][2]
        new_geo::Geodesic = get_geodesic_nocommon_edges(subtree_a, subtree_b)
        geo.ratio_seq = interleave(geo.ratio_seq, new_geo.ratio_seq)
    end # for
    return geo
end # geodesic


"""
    get_common_edge_lengths(trees::Vector{T}, common_edges::Vector{CommonEdge}, 
        l::Int64)::Tuple{Vector{Float64}, Vector{Float64}} where T<:GeneralNode
"""
function get_common_edge_lengths(trees::Vector{T}, common_edges::Vector{CommonEdge}, 
                                l::Int64)::Vector{EdgeLengths} where T<:GeneralNode

    c_e_lengths::Tuple{Vector{Float64}, Vector{Float64}} = ([], [])
    post_orders = post_order.(trees)
    bps = get_bipartitions_as_bitvectors.(trees)
    inc_length::Float64 = 0.0
    ind::Int64 = 0

    for common_edge in common_edges 
        split = get_split(common_edge[1], l)
        for i in [1,2]
            ind = findfirst(x -> x == split, bps[i])
            inc_length = isnothing(ind) ? 0.0 : post_orders[i][ind].inc_length
            println(i)
            push!(c_e_lengths[i], inc_length)
        end # for
    end # for
    return [(c_e_lengths[1][i], c_e_lengths[2][i]) for i in 1:length(c_e_lengths[1])]
end # get_common_edge_lengths


"""
    split_on_common_edge(tree1::T, tree2::T; non_common_edges=[])
        ::Vector{Tuple{T, T}} where T<:GeneralNode

This function finds the common edges of two trees, and splits them at the first common edge
it finds. Then it recursively splits the resulting subtrees aswell. 

Returns a Vector of all pairs of subtrees that share no common edges.

* `tree1` : root node of the first tree

* `tree2` : root node of the second tree

* `non_common_edges` : vector of non common edges. It is initialized empty and appended to 
                       at each recursion
"""
function split_on_common_edge(tree1::T, tree2::T; non_common_edges=[]
                             )::Vector{Tuple{T, T}} where T<:GeneralNode

    trees::Tuple{T, T} = (tree1, tree2)
    num_nodes::Vector{Int64} = [length(post_order(t)) for t in trees]
    leaves::Vector{Vector{T}} = [sort!(get_leaves(t), by=x->x.name) for t in trees]
    num_edges::Vector{Int64} = [num_nodes[i] - length(leaves[i]) - 1 for i in 1:2]
    (num_edges[1] <= 0 || num_edges[2] <= 0) && return []
    
    common_edges::Vector{CommonEdge} = get_common_edges(trees...)
    # if there are no common edges, add the trees to the array of subtrees that share no 
    # common edges
    if isempty(common_edges)
        push!(non_common_edges, trees)
        return non_common_edges
    end # if
    
    # get the first common edge that was found
    common_edge::CommonEdge = common_edges[1]
    # get a bit vector representing the split of the common edge 
    split::BitVector = get_split(common_edge[1], length(leaves[1]))
    
    # find the common node in each tree by using its split
    common_node1, rev1 = get_node_from_split(tree1, split, leaves[1])
    common_node2, rev2 = get_node_from_split(tree2, split, leaves[2])

    # track if we have to swap the tree pairs after using the reversed split for one of 
    # the common nodes 
    reverse::Bool = rev1 ⊻ rev2
    
    # split the trees at the common nodes
    split_tree!(common_node1)
    split_tree!(common_node2)

    # initialize the new subtrees
    initialize_tree!(common_node1)
    initialize_tree!(common_node2)
    initialize_tree!.(trees)
    
    # recursively split new subtrees as well
    if reverse
        split_on_common_edge(tree1, common_node2; non_common_edges)
        split_on_common_edge(tree2, common_node1; non_common_edges)
    else
        split_on_common_edge(common_node1, common_node2; non_common_edges)
        split_on_common_edge(tree1, tree2; non_common_edges)
    end # if/else
    return non_common_edges
end # split_on_common_edge
 

"""
    get_common_edges(tree1::T, tree2::T)::Vector{CommonEdge} where T<:GeneralNode

This function returns all the common edges of two trees with the same leafset. It also
calculates the length difference of each pair of common edges. 

* `tree1` : root node of the first tree

* `tree2` : root node of the second tree
"""
function get_common_edges(tree1::T, tree2::T)::Vector{CommonEdge} where T<:GeneralNode
    common_edges::Vector{CommonEdge} = []
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
            push!(common_edges, (node, length_diff))
        end
    end # for
    return common_edges
end # get_common_edges


"""
    get_node_from_split(tree::T, split::BitVector, leaves::Vector{T}
                       )::Tuple{T, Bool} where T<:GeneralNode

--- INTERNAL ---
This function finds the node in a tree that corresponds to the input split or the reverse of
the input split (in case the found node for the original split would have been the root).

Returns a Vector of all pairs of subtrees that share no common edges.

* `tree` : root node of the tree.

* `split` : A bit vector representing the split

* `leaves` : vector of the leaves of the tree
"""
function get_node_from_split(tree::T, split::BitVector, leaves::Vector{T}=[]
                            )::Tuple{T, Bool} where T<:GeneralNode

    if isempty(leaves)
        leaves = get_leaves(tree)
    end # if
    common_node::T = find_lca(tree, leaves[split])
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
    split_tree!(node::T)::Nothing where T<:GeneralNode

This function splits a tree at the input node. The input node acts as the root of one of the
resulting subtrees. The tree headed by the original root, receives a new node where the
split node used to be. It represents all the leaves that were below the split node.

* `node`: the node where the tree will be split.
"""
function split_tree!(node::T)::Nothing where T<:GeneralNode
    mother = node.mother
    # remove the split node from the tree
    remove_child!(mother, node)
    node.root = true
    # add node representing the leaves below the split node
    add_child!(mother, Node(join(sort!([leaf.name for leaf in get_leaves(node)]), " ")))
    return
end # split_tree!


"""
    interleave(rs1::RatioSequence, rs2::RatioSequence)::RatioSequence

--- INTERNAL ---
Interleaves the ratio sequences rs1 and rs2 after combining them to get the ascending ratio
sequence with the min distance. Returns a new ratio sequence.

* `rs1`: First ratio sequence
* `rs2`: Second ratio sequence
"""
function interleave(rs1::RatioSequence, rs2::RatioSequence)::RatioSequence
    combined1::RatioSequence = get_non_des_rs_with_min_dist(rs1)
    combined2::RatioSequence = get_non_des_rs_with_min_dist(rs2)

    interleaved_rs = RatioSequence()
    ind1::Int64 = 1
    ind2::Int64 = 1
    while (ind1 <= length(combined1.ratios) && ind2 <= length(combined2.ratios))
        if get_ratio(combined1.ratios[ind1]) <= get_ratio(combined2.ratios[ind2])
            push!(interleaved_rs.ratios, combined1.ratios[ind1])
            ind1 += 1
        else
            push!(interleaved_rs.ratios, combined2.ratios[ind2])
            ind2 += 1
        end # if/else
    end # while

    while ind1 <= length(combined1.ratios)
        push!(interleaved_rs.ratios, combined1.ratios[ind1])
        ind1 += 1
    end # while

    while ind2 <= length(combined2.ratios)
        push!(interleaved_rs.ratios, combined2.ratios[ind2])
        ind2 += 1
    end # while

    return interleaved_rs
end # interleave


function get_non_des_rs_with_min_dist(rs::RatioSequence)::RatioSequence
    length(rs.ratios) < 2 && return rs
    combined_rs::RatioSequence = deepcopy(rs)
    i::Int64 = 1
    combine_code::Int64 = 0
    combined_ratio::Ratio = Ratio()
    cc_array::Vector{Int64} = zeros(Int64, length(rs.ratios) - 1) .+ 2
    a::Int64 = 1

    while i < length(combined_rs.ratios) - 1
        if get_ratio(combined_rs.ratios[i]) > get_ratio(combined_rs.ratios[i+1])
            combined_ratio = combine(combined_rs.ratios[i], combined_rs.ratios[i+1])
            deleteat!(combined_rs.ratios,[i, i+1])
            insert!(combined_rs.ratios, i, combined_ratio.ratios)
            cc_array[a] = 1
            if i > 1
                i -= 1
                while cc_array[a] == 1
                    a -= 1
                end # while
            else
                while (a < length(rs.ratios) - 1) && cc_array[a] != 2
                    a += 1
                end # while
            end # if/else
        else
            cc_array[a] = 0
            i += 1
            while (a < length(rs.ratios) - 1) && cc_array[a] != 2
                a += 1
            end # while
        end # if/else
    end # while
    
    for k in 1:(length(rs.ratios) - 1)
        if cc_array[k] == 1
            combine_code = combine_code + 2 ^ k
        end # if
    end # for
    combined_rs.combine_code = combine_code
    return combined_rs
end # get_non_des_rs_with_min_dist


"""
    get_geodesic_nocommon_edges(tree1::T, tree2::T)::Geodesic where T<:GeneralNode
"""
function get_geodesic_nocommon_edges(tree1::T, tree2::T)::Geodesic where T<:GeneralNode
    trees::Tuple{T, T} = (tree1, tree2)
    num_nodes::Vector{Int64} = [length(post_order(t)) for t in trees]
    leaves::Vector{Vector{T}} = [sort!(get_leaves(t), by=x->x.name) for t in trees]
    num_edges::Vector{Int64} = [num_nodes[i] - length(leaves[i]) - 1 for i in 1:2]
    rs::RatioSequence = RatioSequence()
    a_vertices::Vector{Int64} = []
    b_vertices::Vector{Int64} = []
    queue::Vector{Ratio} = Vector{Ratio}()
    ratio::Ratio = Ratio()
    r1::Ratio, r2::Ratio = [Ratio() for _ in 1:2]
    cover::Matrix{Int64} = zeros(Int64, 1, 1)

    common_edges = get_common_edges(trees...)
    # doublecheck to make sure the trees have no common edges
    length(common_edges) != 0 && throw(ArgumentError("Exiting: Can't compute geodesic between subtrees that have common edges."))
    
    if num_edges[1] == 0 || num_edges[2] == 0
        throw(ArgumentError("Exiting: Can't compute the geodesic for trees with no edges"))
    end # if

    if num_edges[1] == 1 || num_edges[2] == 1
        internal_nodes1 = filter!(x -> x.nchild != 0, post_order(tree1)[1:end-1])
        internal_nodes2 = filter!(x -> x.nchild != 0, post_order(tree1)[1:end-1])
		push!(rs.ratios, Ratio(internal_nodes1, internal_nodes2))
		return Geodesic(rs)
	end # if

    leaf_dict::Dict{String, Int64} = Dict(leaf.name => leaf.num for leaf in leaves[1])
    int_nodes1 = filter!(x -> x.nchild != 0, post_order(tree1)[1:end-1])
    int_nodes2 = filter!(x -> x.nchild != 0, post_order(tree2)[1:end-1])
    incidence_matrix::BitMatrix =  get_incidence_matrix(int_nodes1, int_nodes2, leaf_dict)

    graph::BipartiteGraph = build_bipartite_graph(incidence_matrix,
                                                  [n.inc_length for n in int_nodes1], 
                                                  [n.inc_length for n in int_nodes2])
    
    push!(queue, Ratio(int_nodes1, int_nodes2)) 
	
    while length(queue) > 0
        ratio = popfirst!(queue)
        a_vertices = []
        b_vertices = []
        # sizehint!(a_vertices, length(int_nodes1))
        # sizehint!(b_vertices, length(int_nodes2))
        for i in 1:length(ratio.e_edges)
            push!(a_vertices, findfirst(x -> x.num == ratio.e_edges[i].num, int_nodes1)) 
        end # for
        for i in 1:length(ratio.f_edges)
            push!(b_vertices, findfirst(x -> x.num == ratio.f_edges[i].num, int_nodes2)) 
        end # for
        
        cover = get_vertex_cover(graph, a_vertices, b_vertices)
        
        if (cover[1, 1] == 1 || (cover[1, 1] == length(a_vertices))) 
            push!(rs.ratios, ratio)
        
        else
            r1 = Ratio()
            r2 = Ratio()
            j = 1

            for i in 1:length(a_vertices)
                if j < (findfirst(cover[3,:] .== 0) - 1) && (a_vertices[i] == cover[3, j])
                    add_e_edge!(r1, int_nodes1[a_vertices[i]])
                    j += 1
                else
                    add_e_edge!(r2, int_nodes1[a_vertices[i]])
                end # if/else
            end # for
            j = 1

            for i in 1:length(b_vertices)
                if j < (findfirst(cover[4,:] .== 0) - 1) && (b_vertices[i] == cover[4, j])
                    add_f_edge!(r2, int_nodes2[b_vertices[i]])
                    j += 1
                else
                    add_f_edge!(r1, int_nodes2[b_vertices[i]])
                end # if/else
            end # for
        pushfirst!(queue, r2)
        pushfirst!(queue, r1)
        end # if/else
    end # while 
    return Geodesic(rs)     
end # get_geodesic_nocommon_edges


function get_vertex_cover(bg::BipartiteGraph, a_index::Vector{Int64}, b_index::Vector{Int64}
                     )::Matrix{Int64}

    n_AVC = length(a_index)
    n_BVC = length(b_index)
    total::Float64 = 0.0
    ab_flow::Matrix{Float64} = zeros(bg.nums[1], bg.nums[2])
    k, a_scanlistsize, b_scanlistsize, a_pathnode, b_pathnode = [1 for _ = 1:5]
    augmenting_pathend::Int64 = -1
    cd::Matrix{Int64} = zeros(Int64, 4, bg.nums[3])   
    a_scanlist::Vector{Int64} = [0 for _ = 1:bg.nums[1]]
    b_scanlist::Vector{Int64} = [0 for _ = 1:bg.nums[2]]

    # set normalized weights
    for i in 1:n_AVC total += bg.a_vertices[a_index[i]].weight end
    for i in 1:n_AVC
        bg.a_vertices[a_index[i]].residual = bg.a_vertices[a_index[i]].weight / total
    end # for
    total = 0.0
    for i in 1:n_BVC total += bg.b_vertices[b_index[i]].weight end
    for i in 1:n_BVC
        bg.b_vertices[b_index[i]].residual = bg.b_vertices[b_index[i]].weight / total
    end # for

    ### FLOW ALGORITHM ###
    total = 1.0
    while(total > 0.0)
        # scan phase
        # set labels
        total = 0.0
        for i in 1:n_AVC
            bg.a_vertices[a_index[i]].label = -1.0
            bg.a_vertices[a_index[i]].pred = -1.0
        end # for
        for i in 1:n_BVC
            bg.b_vertices[b_index[i]].label = -1.0
            bg.b_vertices[b_index[i]].pred = -1.0
        end # for
        a_scanlistsize = 1
        for i in 1:n_AVC
            if bg.a_vertices[a_index[i]].residual > 0.0
                bg.a_vertices[a_index[i]].label = bg.a_vertices[a_index[i]].residual
                a_scanlist[a_scanlistsize] = a_index[i]
                a_scanlistsize += 1 
            else
                bg.a_vertices[a_index[i]].label = -1.0
            end # if/else
        end # for
        for i in 1:n_BVC bg.b_vertices[i].label = -1.0 end

        # scan for an augmenting path
        a_scanlistsize -= 1
        while(a_scanlistsize != 1)
            # scan the a side nodes
            b_scanlistsize = 1
            for i in 1:(a_scanlistsize)
                for j in 1:n_BVC
                    if bg.incidence_matrix[a_scanlist[i], b_index[j]] && bg.b_vertices[b_index[j]].label == -1.0
                        bg.b_vertices[b_index[j]].label = bg.a_vertices[a_scanlist[i]].label
                        bg.b_vertices[b_index[j]].pred = a_scanlist[i]
                        b_scanlist[b_scanlistsize] = b_index[j]
                        b_scanlistsize += 1
                    end # if
                end # for j
            end # for i
            
            b_scanlistsize -= 1
            # scan the b side nodes
            a_scanlistsize = 1
            for j in 1:b_scanlistsize
                if bg.b_vertices[b_scanlist[j]].residual > 0.0
                    total = min(bg.b_vertices[b_scanlist[j]].residual, 
                                bg.b_vertices[b_scanlist[j]].label)
                    augmenting_pathend = b_scanlist[j]
                    @goto escape_inner_while_loop
                else
                    for i in 1:n_AVC
                        if (bg.incidence_matrix[a_index[i],b_scanlist[j]] && 
                            bg.a_vertices[a_index[i]].label == -1 &&
                            ab_flow[a_index[i], b_scanlist[j]] > 0)
                            bg.a_vertices[a_index[i]].label = 
                                min(bg.b_vertices[b_scanlist[j]].label, 
                                    ab_flow[a_index[i], b_scanlist[j]])

                            bg.a_vertices[a_index[i]].pred = b_scanlist[j]
                            a_scanlist[a_scanlistsize] = a_index[i]
                            a_scanlistsize += 1
                        end # if
                    end # for
                end # if/else
            end # for
        end # while
        @label escape_inner_while_loop
        
        # flow augmentation
        if total > 0.0
            bg.b_vertices[augmenting_pathend].residual -= total
            b_pathnode = augmenting_pathend
            a_pathnode = bg.b_vertices[b_pathnode].pred
            ab_flow[a_pathnode, b_pathnode] += total
            while bg.a_vertices[a_pathnode].pred != -1
                b_pathnode = bg.a_vertices[a_pathnode].pred
                ab_flow[a_pathnode, b_pathnode] -= total
                a_pathnode = bg.b_vertices[b_pathnode].pred
                ab_flow[a_pathnode, b_pathnode] += total
            end # while
            bg.a_vertices[a_pathnode].residual -= total
         else
            k = 1
            for i in 1:n_AVC
                if bg.a_vertices[a_index[i]].label == -1
                    cd[3, k] = a_index[i]
                    k += 1
                end # if
            end # for

            cd[1, 1] = k
            k = 1
            for j in 1:n_BVC
                if bg.b_vertices[b_index[j]].label >= 0
                    cd[4, k]=b_index[j]
                    k += 1
                end # if
            end # for
            
            cd[2, 1] = k
         end # if/else    
    end # while
    return cd
end # get_vertex_cover

"""
    get_incidence_matrix(edges1::Vector{T}, edges2::Vector{T},
        leaf_dict::Dict{String,Int64})::BitMatrix where T<:GeneralNode

--- INTERNAL ---
This function computes the incidence_matrix of two trees. A field of the matrix contains a 1
if the corresponding subtrees of our input trees are not compatible.

Returns the incidence_matrix.

* `edges1` : edges of the first tree

* `edges2` : edges of the second tree

* `leaf_dict` : Dictionary that links each leaf to a unique number
"""
function get_incidence_matrix(edges1::Vector{T}, edges2::Vector{T}, 
                              leaf_dict::Dict{String,Int64})::BitMatrix where T<:GeneralNode

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

* `b1` : first split

* `b2` : second split
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


"""
    geo_avg(edges::Vector{T})::Float64 where T<:GeneralNode

--- INTERNAL ---
This function sums the squared edge lengths of the input edges and returns the square root
of that sum.

* `edges` : root node of the first tree
"""
function geo_avg(edges::Vector{T})::Float64 where T<:GeneralNode
    avg::Float64 = 0.0
    for edge in edges
        avg += edge.inc_length ^ 2
    end # for
    return sqrt(avg)
end # geo_avg


"""
    add_e_edge!(ratio::Ratio, node::T) where T<:GeneralNode

--- INTERNAL ---
This function adds a new edge to one side of a ratio and adjusts the ratio length 
accordingly.

* `ratio` : ratio to which node is added

* `node` : node that is added to the ratio
"""
function add_e_edge!(ratio::Ratio, node::T) where T<:GeneralNode
    push!(ratio.e_edges, node)
    ratio.e_length = sqrt(ratio.e_length ^ 2 + node.inc_length ^ 2)
end # add_e_edge!


"""
    add_f_edge!(ratio::Ratio, node::T) where T<:GeneralNode

--- INTERNAL ---
Works like add_e_edge!
"""
function add_f_edge!(ratio::Ratio, node::T) where T<:GeneralNode
    push!(ratio.f_edges, node)
    ratio.f_length = sqrt(ratio.f_length ^ 2 + node.inc_length ^ 2)
end # add_f_edge!


"""
    get_ratio(r::Ratio)::Float64

--- INTERNAL ---
Computes the ratio of a ratio struct

* `ratio` : Ratio for which ratio is computed
"""
function get_ratio(r::Ratio)::Float64
    r.e_length / r.f_length
end # get_ratio


function combine(r1::Ratio, r2::Ratio)::Ratio
    r::Ratio = Ratio()
    edges::Vector{T where T<:GeneralNode} = []
    if length(r1.e_edges) == 0 && length(r2.e_edges) == 0
        r.e_length(sqrt(r1.e_length ^ 2 + r2.e_length ^ 2))
    else
        edges = deepcopy.([r1.e_edges, r2.e_edges])
        addall_e_edges!.(r, edges)
    end # if/else

    if length(r1.f_edges) == 0 && length(r2.f_edges) == 0
        r.f_length(sqrt(r1.f_length ^ 2 + r2.f_length ^ 2))
    else 
        edges = deepcopy.([r1.f_edges, r2.f_edges])
        addall_f_edges!.(r, edges)
    end # if/else
    return r
end # combine

function addall_e_edges!(r::Ratio, edges::Vector{T}) where T<:GeneralNode
    r.e_edges = edges
    r.e_length = geo_avg(r.e_edges) 
end # addall_e_edges!

function addall_f_edges!(r::Ratio, edges::Vector{T}) where T<:GeneralNode
    r.f_edges = edges
    r.f_length = geo_avg(r.f_edges) 
end # addall_f_edges!

