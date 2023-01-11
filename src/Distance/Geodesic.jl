"""
    Geodesic
    
Struct that tracks all relevant parameters, when computing the Geodesic between two trees.

* `ratio_seq` : The ratio sequence of the Geodesic.
* `leaf_contribution` : Contribution of the leafs to the geodesic distance.
* `common_edges` : Stores the common edges of the two compared trees.
"""
mutable struct Geodesic{T<:GeneralNode}
    ratio_seq::RatioSequence
    leaf_contribution²::Float64
    common_edges::Vector{Tuple{T, T}}

    Geodesic(rs::RatioSequence) = 
        new{GeneralNode}(rs, 0.0, [])
end # Geodesic


"""
    get_distance(geo::Geodesic)::Float64

--- INTERNAL ---
Returns the distance of a geodesic. The smaller the number, the more similar are the trees.

* `geo` : The geodesic for which the distance is calculated.
"""
function get_distance(geo::Geodesic, verbose::Bool)::Float64
    common_edge_dist²::Float64 = common_edge_contribution(geo.common_edges)
    verbose && println("Common edge contribution squared: $common_edge_dist²")
    non_descending_rs::RatioSequence = get_non_desc_rs_with_min_dist(geo.ratio_seq)
    return sqrt(get_distance(non_descending_rs) ^ 2 + common_edge_dist² + geo.leaf_contribution²)
end # get_distance


"""
    common_edge_contribution(common_edges::Vector{Tuple{T}})::Float64 where T<:GeneralNode

--- INTERNAL ---
Returns the squared common edge contribution of a vector of common edges.

* `common_edges` : Common edges of the two trees.
"""
function common_edge_contribution(common_edges::Vector{Tuple{T, T}})::Float64 where T<:GeneralNode
    common_edge_dist²::Float64 = 0.0
    for ce_pair in common_edges
        common_edge_dist² += (ce_pair[1].inc_length - ce_pair[2].inc_length) ^ 2
    end # for
    return common_edge_dist²
end # common_edge_contribution


"""
    Vertex
    
This struct represents a vertex of a bipartite graph.
"""
mutable struct Vertex
    label::Float64
    weight::Float64
    residual::Float64
    pred::Int64

    Vertex(weight::Float64) = new(0.0, weight*weight, 0.0, 0)
end # Vertex


"""
    BipartiteGraph
    
This struct represents a bipartite node-weighted graph.

* `incidence_matrix` : The incidence_matrix of the graph.
* `nums` : The number of nodes on each side, and their maximum.
* `a_vertices` : A vector of the vertices of the first tree.
* `b_vertices` : A vector of the vertices of the second tree.
"""
mutable struct BipartiteGraph
    incidence_matrix::BitMatrix
	nums::Vector{Int64}
	a_vertices::Vector{Vertex}
    b_vertices::Vector{Vertex}
end # BipartiteGraph


"""
    build_bipartite_graph(incidence_matrix::BitMatrix, a_weight::Vector{Float64}, 
                          b_weight::Vector{Float64})::BipartiteGraph
  
--- INTERNAL ---                          
This function builds and returns bipartite graph based on an incidence matrix and node
weights given to it.

`incidence_matrix` : The incidence matrix of the graph.
`a_weights` : Weights of the a-side nodes.
`b_weights` : Weights of the b-side nodes.
"""
function build_bipartite_graph(incidence_matrix::BitMatrix, a_weights::Vector{Float64}, 
                               b_weights::Vector{Float64})::BipartiteGraph

    a::Int64 = length(a_weights)
    b::Int64 = length(b_weights)
    n::Int64 = max(a, b)
    a_vertices::Vector{Vertex} = fill(Vertex(0.0), n)
    b_vertices::Vector{Vertex} = fill(Vertex(0.0), n)
    for i in 1:a
        a_vertices[i]= Vertex(a_weights[i])
    end # for 
    for i in 1:b
        b_vertices[i]= Vertex(b_weights[i])
    end # for
    return BipartiteGraph(incidence_matrix, [a, b, n], a_vertices, b_vertices)
end # build_bipartite_graph


"""
    get_common_edges(tree1::T, tree2::T)
        ::Vector{Tuple{T, T}} where T<:GeneralNode

This function finds the common edges of two trees with the same leafset.

Returns the common edges as a vector of tuples.   
    
* `tree1` : root node of the first tree
* `tree2` : root node of the second tree
"""
function get_common_edges(tree1::T, tree2::T)::Vector{Tuple{T, T}} where T<:GeneralNode

    common_edges::Vector{Tuple{T, T}} = []
    tree_splits::Vector{BitVector} = get_bipartitions_as_bitvectors(tree2)
    # ToDo: Ugly
    l = length(collect(get_leaves(tree1)))
    nodes_tree2::Vector{T} = collect(post_order(tree2))
    sort!(nodes_tree2, by = x -> x.num)
    #bp::Vector{BitVector} = get_bipartitions_as_bitvectors(tree2)

    for node in post_order(tree1)
        (node.nchild == 0 || isroot(node)) && continue
        # get the split of the current node represented as a BitVector
        split::BitVector = get_split(node, l)
        # if the same split exits in both trees, then we found a common node
        for (ind, s_split) in enumerate(tree_splits)
            # find the node in the other tree and save its length
            #ind = findfirst(x -> x == split, bp)
            if s_split == split
                # push the common node of each tree and their lengths
                push!(common_edges, (node, nodes_tree2[ind]))
            end # if
        end 
    end # for
    return common_edges
end # get_common_edges


"""
    split_on_common_edge(tree1::T, tree2::T; non_common_edges=[])
        ::Vector{Tuple{T, T}} where T<:GeneralNode

--- INTERNAL ---
This function finds the common edges of two trees, and splits them at the first common edge
it finds. Then it recursively splits the resulting subtrees aswell. 

Returns a Vector of all pairs of subtrees that share no common edges.

* `tree1` : root node of the first tree.
* `tree2` : root node of the second tree.
* `non_common_edges` : vector of non common edges. It is initialized empty and appended to 
                       at each recursion.
"""
function split_on_common_edge(tree1::T, tree2::T; non_common_edges=[]
                             )::Vector{Tuple{T, T}} where T<:GeneralNode

    trees::Tuple{T, T} = (tree1, tree2)
    num_nodes::Vector{Int64} = [treesize(t) for t in trees]
    #ToDo: Ugly
    leaves::Vector{Vector{T}} = [sort!(collect(get_leaves(t)), by=x->x.name) for t in trees]
    num_edges::Vector{Int64} = [num_nodes[i] - length(leaves[i]) - 1 for i in 1:2]
    (num_edges[1] <= 0 || num_edges[2] <= 0) && return []
    
    common_edges::Vector{Tuple{T, T}} = get_common_edges(trees...)
    # if there are no common edges, add trees to list of trees that share no common edges
    if isempty(common_edges)
        push!(non_common_edges, trees)
        return non_common_edges
    end # if
    
    # get the common nodes
    common_node1 = common_edges[1][1]
    common_node2 = common_edges[1][2]

    # split the trees at the common nodes
    split_tree!(common_node1)
    split_tree!(common_node2)

    # initialize the new subtrees
    initialize_tree!(common_node1)
    initialize_tree!(common_node2)
    initialize_tree!.(trees)
    
    # recursively split new subtrees as well
    split_on_common_edge(common_node1, common_node2; non_common_edges)
    split_on_common_edge(tree1, tree2; non_common_edges)
    return non_common_edges
end # split_on_common_edge


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
    get_geodesic_nocommon_edges(tree1::T, tree2::T)::Geodesic where T<:GeneralNode

--- INTERNAL ---
This function computes the geodesic between two trees that share no common edges. It is 
called in the function `geodesic` after the common edges have been removed.
"""
function get_geodesic_nocommon_edges(tree1::T, tree2::T)::Geodesic where T<:GeneralNode
    trees::Tuple{T, T} = (tree1, tree2)
    num_nodes::Vector{Int64} = [treesize(t) for t in trees]
    leaves::Vector{Vector{T}} = [sort!(collect(get_leaves(t)), by=x->x.name) for t in trees]
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
        # ToDo: Ugly
        internal_nodes1 = filter!(x -> x.nchild != 0, collect(post_order(tree1))[1:end-1])
        internal_nodes2 = filter!(x -> x.nchild != 0, collect(post_order(tree2))[1:end-1])
		push!(rs.ratios, Ratio(internal_nodes1, internal_nodes2))
		return Geodesic(rs)
	end # if

    leaf_dict::Dict{String, Int64} = Dict(leaf.name => leaf.num for leaf in leaves[1])
    # ToDo: Ugly
    int_nodes1 = filter!(x -> x.nchild != 0, collect(post_order(tree1))[1:end-1])
    int_nodes2 = filter!(x -> x.nchild != 0, collect(post_order(tree2))[1:end-1])
    incidence_matrix::BitMatrix =  get_incidence_matrix(int_nodes1, int_nodes2, leaf_dict)

    graph::BipartiteGraph = build_bipartite_graph(incidence_matrix,
                                                  [n.inc_length for n in int_nodes1], 
                                                  [n.inc_length for n in int_nodes2])
    
    push!(queue, Ratio(int_nodes1, int_nodes2)) 
	
    while length(queue) > 0
        ratio = popfirst!(queue)
        a_vertices = []
        b_vertices = []
        for i in 1:length(ratio.e_edges)
            push!(a_vertices, findfirst(x -> x.num == ratio.e_edges[i].num, int_nodes1)) 
        end # for
        for i in 1:length(ratio.f_edges)
            push!(b_vertices, findfirst(x -> x.num == ratio.f_edges[i].num, int_nodes2)) 
        end # for
        
        cover = get_vertex_cover(graph, a_vertices, b_vertices)
        
        if (cover[1, 1] == 1 || (cover[1, 1] == length(a_vertices) + 1)) 
            push!(rs.ratios, ratio)
        
        else
            r1 = Ratio()
            r2 = Ratio()
            j = 1

            for i in 1:length(a_vertices)
                if j < length(cover[3,:]) && (a_vertices[i] == cover[3, j])
                    add_e_edge!(r1, int_nodes1[a_vertices[i]])
                    j += 1
                else
                    add_e_edge!(r2, int_nodes1[a_vertices[i]])
                end # if/else
            end # for
            j = 1

            for i in 1:length(b_vertices)
                if j < length(cover[4,:]) && (b_vertices[i] == cover[4, j])
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


"""
    get_vertex_cover(bg::BipartiteGraph, a_index::Vector{Int64}, b_index::Vector{Int64})
        ::Matrix{Int64}

--- INTERNAL ---
This function computes the min-normalized-square-weighted vertex cover. 

Returns a 4xn matrix cd[1][] = #A-side cover elements
                     cd[2][] = #B-side cover elements
                     cd[3][] = list of A-side cover elements
                     cd[4][] = list of B-side cover elements

* `bg` : bipartite graph, holds incidence_matrix, several number used for calculations, 
         as well as the a_vertices and b_vertices.
* `a_index` : Vector of indices for a_vertices.
* `b_index` : Vector of indices for b_vertices.
"""
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
        while(a_scanlistsize != 1)
            # scan the a side nodes
            b_scanlistsize = 1
            for i in 1:(a_scanlistsize - 1)
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
            for j in 1:(b_scanlistsize - 1)
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
                    cd[4, k] = b_index[j]
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

* `edges1` : edges of the first tree.
* `edges2` : edges of the second tree.
* `leaf_dict` : Dictionary that links each leaf to a unique number.
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
"""
function crosses(b1::BitVector, b2::BitVector)::Bool
    disjoint::Bool = all(.!(b1 .& b2))
    contains1::Bool = all(b1 .>= b2)
    contains2::Bool = all(b2 .>= b1)
    return !(disjoint || contains1 || contains2)   
end # crosses


"""
    geo_avg(edges::Vector{T})::Float64 where T<:GeneralNode

--- INTERNAL ---
This function sums the squared edge lengths of the input edges and returns the square root
of that sum.
"""
function geo_avg(edges::Vector{T})::Float64 where T<:GeneralNode
    avg::Float64 = 0.0
    for edge in edges
        avg += edge.inc_length ^ 2
    end # for
    return sqrt(avg)
end # geo_avg
