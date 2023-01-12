"""
    prune_tree(root::T, node_names::Vector{String})::T where T<:AbstractNode

This function removes specific nodes, including their descendants, from a tree.

* `root` : root Node of tree to prune.

* `node_names` : vector of strings, used to specify nodes to remove.
"""
function prune_tree(root::T, node_names::Vector{String})::T where T<:AbstractNode
    # copy the tree and call the inplace version of the function on the copy
    if root.name in node_names
        throw(ArgumentError("trying to prune root, please set root to nothing instead"))
    end
    copyroot = deepcopy(root)
    prune_tree!(copyroot, node_names)
    return copyroot
end


"""
    prune_tree!(root::T, node_names::Vector{String})::Nothing where T<:AbstractNode

In-place version of prune_tree.

* `root` : root Node of tree to prune.

* `node_names` : vector of strings, used to specify nodes to remove.
"""
function prune_tree!(root::T, node_names::Vector{String})::Nothing where T<:AbstractNode
    nodes = post_order(root)
    names = [node.name for node in nodes]
    nodes_to_prune = Vector{String}()
    for name in node_names
        indices = findall(x->x == name, names)
        if length(indices) > 1
            throw(ArgumentError("Multiple nodes are named \"$name\""))
        elseif length(indices) == 0
            print("Warning: No node named \"$name\" in tree")
        else
            push!(nodes_to_prune, name)
        end # if/else
    end # for
    if length(nodes_to_prune) == 0
        throw(ArgumentError("None of the node names correspond to a node in the tree"))
    end # if
    prune_tree!(root, find_by_name.(Ref(root), nodes_to_prune))
end


"""
    prune_tree!(root::T, node_names::Vector{T})::Nothing where T<:AbstractNode

In-place version of prune_tree.

* `root` : root node of tree to prune.

* `node_names`: vector of Node objects to be removed from tree.
"""
function prune_tree!(root::T, nodes::Vector{T})::Nothing where T<:AbstractNode
    if root in nodes
        throw(ArgumentError("trying to prune root, please set root to nothing instead"))
    else
        for node in nodes
            remove_child!(get_mother(node), node)
        end
    end
end
