
#################### Generics ####################
const fulltreeit = Union{PostOrderDFS, PreOrderDFS, StatelessBFS}

Base.length(n::T) where T<:fulltreeit = treesize(n.root)
Base.lastindex(n::T) where T<:fulltreeit = treesize(n.root)
Base.getindex(n::T, itr) where T<:TreeIterator{<:GeneralNode} = [n[i] for i in itr]
function Base.getindex(n::T, ind::Int) where T<:TreeIterator{<:GeneralNode}
    ct = 1
    for node in n
        ind == ct && return node
        ct += 1
    end
    throw(BoundsError("Index $ind is out of bounds for tree with $(treesize(n.root)) nodes."))
end

#################### Post order traversal ####################

"""
    post_order(root::T) where T<:GeneralNode

This function does post order traversal. Only the root node needs to be supplied.

Returns a post order iterator.

* `root` : root Node of tree.
"""
function post_order(root::T) where T<:GeneralNode
    return PostOrderDFS(root)
end # function post_order


"""
    get_leaves(root::T) where T<:AbstractNode

This function returns the leaves of a tree. Only the root node needs to be supplied.

Returns an Iterator over leaf Nodes.

* `root` : root Node of tree.
"""
function get_leaves(root::T) where T<:GeneralNode
    Leaves(root)
end # function post_order



#################### Pre order traversal ####################


"""
    pre_order(root::T) where T<:GeneralNode

This function does pre order traversal. Only the root node needs to be supplied.

Returns a preoder Iterator.

* `root` : root Node of tree.
"""
function pre_order(root::T) where T<:GeneralNode
    PreOrderDFS(root)
end # function pre_order

#################### Level order traversal ####################

"""
    level_order(node::T) where T<:GeneralNode

This function does level order traversal. Only the root node needs to be supplied.

Returns a level order iterator.

* `node` : root Node of tree.
"""
function level_order(node::T) where T<:GeneralNode
    StatelessBFS(node)
end # function level_order
