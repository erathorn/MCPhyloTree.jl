mutable struct Ratio
    e_length::Float64
	f_length::Float64
	e_edges::Vector{FNode}
	f_edges::Vector{FNode}

    Ratio() = new(0.0, 0.0, [], [])
    Ratio(e::Vector{FNode}, f::Vector{FNode}) = new(geo_avg(e), geo_avg(f), e, f)
end # Ratio

const RatioSequence = Vector{Ratio}
const EdgeLengths = Tuple{Float64, Float64}
const CommonEdge = Tuple{FNode, Float64}

mutable struct Geodesic
    ratio_seq::RatioSequence
    e_leaf_attribs::Vector{Float64}
    f_leaf_attribs::Vector{Float64}
    leaf_contribution²::Float64
    common_edges::CommonEdge
    common_edge_lengths::EdgeLengths

    Geodesic(rs::RatioSequence, e_lengths::Vector{Float64}, f_lengths::Vector{Float64}) = 
        new(rs, e_lengths, f_lengths, 0.0, CommonEdge[], EdgeLengths[])

    Geodesic(rs::RatioSequence) = 
        new(rs, [], [], 0.0, CommonEdge[], EdgeLengths[])
end # Geodesic


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


function geodesic(tree1::FNode, tree2::FNode)
    leaf_contribution²::Float64 = 0.0
    leaves::Vector{FNode} = get_leaves(tree1)
    leaves2::Vector{FNode} = get_leaves(tree2)
    perm = sortperm([leaf.name for leaf in leaves])
    perm2 = sortperm([leaf.name for leaf in leaves2])
    leaf_vec::Vector{String} = [leaf.name for leaf in leaves][perm]
    leaf_vec2::Vector{String} = [leaf.name for leaf in leaves2][perm2]
    leaf_attribs::Vector{Float64} = [leaf.inc_length for leaf in leaves][perm]
    leaf_attribs2::Vector{Float64} = [leaf.inc_length for leaf in leaves2][perm2]
    
    leaf_vec != leaf_vec2 && 
        throw(ArgumentError("The two input trees do not have the same sets of leaves")) 
    for (leaf1, leaf2) in zip(leaves1, leaves2)
        leaf_contribution² += abs(leaf1.inc_length - leaf2.inc_length) ^ 2
    end # for
    
    geo = Geodesic(RatioSequence(), leaf_attribs, leaf_attribs2)
    geo.leaf_contribution² = leaf_contribution²

    non_common_edges::Vector{Tuple{FNode, FNode}} = split_on_commonedge(deepycopy(tree1), 
                                                                      deepcopy(tree2))

    common_edges::Vector{CommonEdge} = get_commonedges(tree1, tree2)
    geo.common_edges = common_edges

    c_e_lengths::Vector{EdgeLengths} = get_commonedge_lengths([tree1, tree2], common_edges,
                                                              length(leaves1))                                                     
    geo.common_edge_lengths = c_e_lengths

    for i in 1:length(non_common_edges)
        subtree_a = non_common_edges[i][1]
        subtree_b = non_common_edges[i][2]
        new_geo::Geodesic = get_geodesic_nocommonedges(subtree_a, subtree_b)
        #TODO: interleave function
        geo.ratio_seq = interleave(geo.ratio_seq, new_geo.ratio_seq)
    end # for
    return geo
end # geodesic


"""
    get_commonedge_lengths(trees::Vector{FNode}, common_edges::Vector{CommonEdge}, 
        l::Int64)::Tuple{Vector{Float64}, Vector{Float64}}
"""
function get_commonedge_lengths(trees::Vector{FNode}, common_edges::Vector{CommonEdge}, 
                                l::Int64)::Vector{EdgeLengths}

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
    lengths::Vector{EdgeLengths} = 
        [(c_e_lengths[1][i], c_e_lengths[2][i]) for i in 1:length(c_e_lengths[1])]
end # get_common_edge_lengths


"""
    split_on_commonedge(tree1::FNode, tree2::FNode; non_common_edges=[])
        ::Vector{Tuple{FNode, FNode}}

This function finds the common edges of two trees, and splits them at the first common edge
it finds. Then it recursively splits the resulting subtrees aswell. 

Returns a Vector of all pairs of subtrees that share no common edges.

* `tree1` : root node of the first tree

* `tree2` : root node of the second tree

* `non_common_edges` : vector of non common edges. It is initialized empty and appended to 
                       at each recursion
"""
function split_on_common_edge(tree1::FNode, tree2::FNode; non_common_edges=[]
                          )::Vector{Tuple{FNode, FNode}}

    numNodes1::Int64 = length(post_order(tree1))
    numNodes2::Int64 = length(post_order(tree2))
    leaves1::Vector{FNode} = sort!(get_leaves(tree1), by=x->x.name)
    leaves2::Vector{FNode} = sort!(get_leaves(tree2), by=x->x.name)
    numEdges1::Int64 = numNodes1 - length(leaves1) - 1
    numEdges2::Int64 = numNodes2 - length(leaves2) - 1
    (numEdges1 <= 0 || numEdges2 <= 0) && return []
    
    common_edges::Vector{CommonEdge} = get_commonedges(tree1, tree2)
    
    # if there are no common edges, add the trees to the array of subtrees that share no 
    # common edges
    if isempty(common_edges)
        push!(non_common_edges, (tree1, tree2))
        return non_common_edges
    end # if
    
    # get the first common edge that was found
    common_edge::CommonEdge = common_edges[1]
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
        split_on_commonedge(tree1, common_node2; non_common_edges)
        split_on_commonedge(tree2, common_node1; non_common_edges)
    else
        split_on_commonedge(common_node1, common_node2; non_common_edges)
        split_on_commonedge(tree1, tree2; non_common_edges)
    end # if/else
    return non_common_edges
end # split_on_commonedge
 

"""
    get_commonedges(tree1::FNode, tree2::FNode)::Vector{CommonEdge}

This function returns all the common edges of two trees with the same leafset. It also
calculates the length difference of each pair of common edges. 

* `tree1` : root node of the first tree

* `tree2` : root node of the second tree
"""
function get_commonedges(tree1::FNode, tree2::FNode)::Vector{CommonEdge}
    commonEdges::Vector{CommonEdge} = []
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
end # get_commonedges


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
    a_vertices::Vector{Int64} = []
    b_vertices::Vector{Int64} = []
    queue::Vector{Ratio} = Vector{Ratio}()
    ratio::Ratio = Ratio()
    r1::Ratio, r2::Ratio = [Ratio() for _ in 1:2]
    cover::Array{Int64} =[[]]

    commonedges = get_commonedges(tree1, tree2)
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
            a_vertices[i] = findfirst(isequal(ratio.e_edges[i]), int_nodes1) 
        end # for
        for i in 1:length(ratio.f_edges)
            b_vertices[i] = findfirst(isequal(ratio.f_edges[i]), int_nodes2) 
        end # for
        
        cover = vertex_cover(bg, a_vertices, b_vertices)
        
        if (cover[1, 1] == 0 || (cover[1, 1] == length(a_vertices))) 
            push(rs, ratio)
        
        else
            r1 = Ratio()
            r2 = Ratio()
            j = 1

            for i in 1:length(a_vertices)
                if j < size(cover)[2] && (a_vertices[i] == cover[3, j])
                    add_e_edge!(r1, find_num(tree1, a_vertices[i]))
                    j += 1
                else
                    add_e_edge!(r2, find_num(tree1, a_vertices[i]))
                end # if/else
            end # for
            j = 1

            for i in 1:b_vertices
                if j < size(cover)[2] && (b_vertices[i] == cover[4, j])
                    add_f_edge!(r2, find_num(tree2, b_vertices[i]))
                    j += 1
                else
                    add_f_edge!(r1, find_num(tree2, b_vertices[i]))
                end # if/else
            end # for
        end # if/else
        pushfirst!(queue, r2)
        pushfirst!(queue, r1)
    end # while 
    return Geodesic(rs)     
end # get_geodesic_nocommonedges


function vertex_cover(bg::BipartiteGraph, a_index::Vector{Int64}, b_index::Vector{Int64}
                     )::Matrix{Int64}

    n_AVC = length(a_index)
    n_BVC = length(b_index)
    total::Float64 = 0.0
    ab_flow::Matrix{Float64} = zeros(n_AVC, n_BVC)
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
        for i in 1:n_BVC b_vertices[i].label = -1.0 end

        #scan for an augmenting path
        while(a_scanlistsize != 1)
            # scan the a side nodes
            b_scanlistsize = 1
            for i in 1:a_scanlistsize
                for j in 1:n_BVC
                    if bg.incidence_matrix[a_scanlist[i], b_index[j]] && bg.b_vertices[b_index[j]].label == -1.0
                        bg.b_vertices[b_index[j]].label = bg.a_vertices[a_scanlist[i]].label
                        bg.b_vertices[b_index[j]].pred = a_scanlist[i]
                        b_scanlist[b_scanlistsize] = b_index[j]
                        b_scanlistsize += 1
                    end # if
                end # for j
            end # for i
            
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
            a_pathnode = bg.b_vertices[b_pathnode].pre_order
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
end # vertex_cover

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


"""
    geo_avg(edges::Vector{FNode})::Float64

--- INTERNAL ---
This function sums the squared edge lengths of the input edges and returns the square root
of that sum.

* `edges` : root node of the first tree
"""
function geo_avg(edges::Vector{FNode})::Float64
    avg::Float64 = 0.0
    for edge in edges
        avg += edge.inc_length ^ 2
    end # for
    return sqrt(avg)
end # geo_avg


"""
    add_e_edge!(ratio::Ratio, node::FNode)

--- INTERNAL ---
This function adds a new edge to one side of a ratio and adjusts the ratio length 
accordingly.

* `ratio` : ratio to which node is added
* `node`: node that is added to the ratio
"""
function add_e_edge!(ratio::Ratio, node::FNode)
    push!(ratio.e_edges, node)
    ratio.e_length = sqrt(e_length ^ 2 + node.inc_length ^ 2)
end # add_e_edge!


"""
    add_f_edge!(ratio::Ratio, node::FNode)

--- INTERNAL ---
Works like add_e_edge!
"""
function add_f_edge!(ratio::Ratio, node::FNode)
    push!(ratio.f_edges, node)
    ratio.f_length = sqrt(f_length ^ 2 + node.inc_length ^ 2)
end # add_f_edge!