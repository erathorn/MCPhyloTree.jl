
"""
to_df(root::GeneralNode)::Tuple{Array{Float64}, Vector{String}}

This function returns a matrix representation of the tree structure and a vector with the column names.
The entry `mat[i,j]` is the length of the edge connecting node `i` with node `j`.
Returns Tuple containing the matrix and a vector of names.

* `root` : root of tree used to create matrix represenation.
"""
function to_df(root::GeneralNode)::Tuple{Array{Float64},Vector{String}}

    #post_order_iteration = post_order(root)

    name_list = [i.name for i in  post_order(root)]
    num_list = [i.num for i in  post_order(root)]
    temp_ar = zeros(Float64, (treesize(root), treesize(root)))
    for i in  post_order(root)
        if i.nchild != 0
            ind = findfirst(isequal(i.num), num_list)
            for j in i.children
                ind2 = findfirst(isequal(j.num), num_list)
                temp_ar[ind[1], ind2[1]] = j.inc_length
            end # end for
        end # end if
    end # end for

    return temp_ar, name_list
end # end function to_df



#################### Covariance wrapper ####################

function to_covariance(tree::N) where {N<:AbstractNode}
    blv = get_branchlength_vector(tree)
    to_covariance(tree, blv)
end # end to_covariance

#################### To Matrix convertes ####################
"""
    to_distance_matrix(tree::T)::Array{Float64,2} where T <:AbstractNode

Calculate the distance matrix over the set of leaves.

Returns an Array of Floats.

* `tree` : root node of tree used to perform caclulcation.
"""
function to_distance_matrix(tree::T)::Array{Float64,2} where {T<:AbstractNode}
    leaves = get_leaves(tree)
    ll = mapreduce(x->1, +, leaves)
    distance_mat = zeros(Float64, ll, ll)
    for i = 1:ll, j = 1:ll
        if i > j
            d = node_distance(tree, leaves[i], leaves[j])
            distance_mat[i, j] = d
            distance_mat[j, i] = d
        end # if
    end #for
    distance_mat
end # function to_distance_matrix

"""
    to_covariance(tree::N, blv::Array{T})::Array{T,2} where {N<:AbstractNode,T<: Real}

Calcualte the variance-covariance matrix from `tree`. An entry (i,j) of the matrix
is defined as the length of the path connecting the latest common ancestor
of i and j with the root of the tree.

Returns an Array of Real numbers.

* `tree` : Node in tree of interest.

* `blv` : branchlength vector of tree.

"""
function to_covariance(tree::N, blv::Vector{T})::Array{T,2} where {N<:AbstractNode,T<:Real}
    leaves::Vector{N} = sort(collect(get_leaves(tree)), by = x->x.num)
    ll = length(leaves)
    covmat = zeros(T, ll, ll)
    
    @inbounds for ind = 1:ll, jnd = 1:ind
        itm = leaves[ind]
        if ind == jnd
            covmat[ind,jnd] =  reduce(+, @view blv[get_path(tree, itm)])
        else
            lca = find_lca(tree, itm, leaves[jnd])
            if !isroot(lca)
                tmp = reduce(+, @view blv[get_path(tree, lca)])
                covmat[ind,jnd] = covmat[jnd,ind] = tmp    
            end # if
        end # if

    end # for
    
    covmat
end# function to_covariance

"""
    function leave_incidence_matrix(root::G)::Matrix{Float64} where {G<:AbstractNode}

Calculate the incidence matrix of the tree whos root node is `root`
For a tree with ``m`` leaves and ``n`` vertecies this function returns an ``m \\times n`` matrix ``L``,
where ``L_{ij} = 1`` if vertex ``j`` is on the path from leave ``i`` to the root of the tree and ``0`` otherwise.

Returns leave incidence matrix.

* `root` : Root node of the tree
"""
function leave_incidence_matrix(root::G)::Matrix{Float64} where {G<:AbstractNode}
    n = treesize(root)-1
    leave_incidence_matrix(root, n)
end


function leave_incidence_matrix(root::G, n::Int)::Matrix{Float64} where {G<:AbstractNode}
    
    out = zeros(maximum(i->i.num, get_leaves(root)), n)
    for leave in get_leaves(root)
        mother = leave
        ln = leave.num
        while !isroot(mother)
            out[ln, mother.num] = 1.0
            mother = get_mother(mother)
        end
    end
    out
end