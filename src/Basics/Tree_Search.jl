
#################### Safe searches ####################

"""
    find_by_binary(tree::GeneralNode, bin::String)

Alias for `find_binary`.
"""
find_by_binary(tree::GeneralNode, bin::String) = find_binary(tree, bin)

"""
    find_by_name(tree::GeneralNode, name::AbstractString)

Alias for `find_name`.
"""
find_by_name(tree::GeneralNode, name::AbstractString) = find_name(tree, name)

"""
    find_by_num(tree::GeneralNode, num::Int64)

Alias for `find_num`.
"""
find_by_num(tree::GeneralNode, num::Int64) = find_num(tree, num)




#################### Unsafe searches ####################

"""
    find_num(root::T, num::Int64)  where T<:GeneralNode

Find a node by its number. The function assumes that the node is present in the
tree.

Do not use this function if you are unsure whether the node is in the tree at all.

Returns reference to Node.

* `root` : root Node of tree to be searched.

* `num` : number of desired Node.
"""
function find_num(root::T, num::I)::T  where {T<:GeneralNode, I<:Integer}
    store = T[]
    found = find_num(root, num, store)
    if length(store) == 0
        throw(ArgumentError("Node not found"))
    else
        return store[1]
    end
end


"""
    find_num(root::T, num::Int64, rn::Vector{T})::Bool  where T<:GeneralNode

Do a post order traversal to find the node corresponding to the `num`.

Returns true if node is found, false otherwise. Desired Node is pushed to rn.

* `root` : root Node of tree to be searched.

* `num` : number of desired Node.

* `rn` : Vector of Nodes; desired Node is pushed to this vector when found.
"""
function find_num(root::T, num::I, rn::Vector{T})::Bool  where {T<:GeneralNode, I<:Integer}
    # if the current node is the correct one store it in rn
    if root.num === num
        push!(rn, root)
        found = true
    else
        found = false
    end

    if !found
        # if the node is not yet found continue
        for child in root.children
            found = find_num(child, num,  rn)
            found && break
        end
    end # if
    return found
end


"""
    find_binary(root::T, bin::String)::T where T<:GeneralNode

Find a node by its binary representation. The function assumes that the node is
present in the tree.

Do not use this function if you are unsure whether the node is in the tree at all.

Returns a reference to the desired Node.

* `root` : root Node of tree to search.

* `bin` : binary representation of desired Node as a String.
"""
function find_binary(root::T, bin::String)::T where T<:GeneralNode
    rv = root
    for i in split(bin[1:end],",")[2:end]
        rv = rv.children[parse(Int64, i)+1]
    end
    rv
end


"""
    find_name(root::T, name::S)::T  where {T<:GeneralNode, S<:AbstractString}

Find a node by its name. Returns reference to Node.

* `root` : root Node of tree to be searched.

* `name` : name of desired Node.
"""
function find_name(root::T, name::S)::T  where {T<:GeneralNode, S<:AbstractString}
    store = T[]
    found = find_name(root, name, store)
    if length(store) == 0
        throw(ArgumentError("Node not found"))
    else
        return store[1]
    end
end

"""
    find_name(root::T, name::S, rn::Vector{T})::Bool where {T<:GeneralNode, S<:AbstractString}

Do a post order traversal to find the node corresponding to the `name`.

Returns true if node is found, false otherwise. Desired Node is pushed to rn.

* `root` : root Node of tree to be searched.

* `name` : name of desired Node.

* `rn` : Vector of Nodes; desired Node is pushed to this vector when found.
"""
function find_name(root::T, name::S, rn::Vector{T})::Bool where {T<:GeneralNode, S<:AbstractString}
    # if the current node is the correct one store it in rn
    if root.name === name
        push!(rn, root)
        found = true
    else
        found = false
    end

    if !found
        # if the node is not yet found continue
        for child in root.children
            found = find_name(child, name,  rn)
            found && break
        end
    end # if
    return found
end


"""
    find_root(node::T)::T where T <: GeneralNode

Finds the root of tree indicated by Node.

Returns reference to root Node of the tree.

* `node` : Node in Tree of interest.
"""
function find_root(node::T)::T where T <: GeneralNode
    while node.root == false
        node = get_mother(node)
    end # while
    return node
end
