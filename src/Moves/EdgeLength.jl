
### Not altering the tree length ###

"""
slide!(root::T) where T<:GeneralNode

This functin performs a slide move on an intermediate node. The node is moved
upwards or downwards on the path specified by its mother and one of its
daughters.

* `root` : root Node of tree.
"""
function slide!(root::T) where {T<:GeneralNode}

    available = [node.num for node in post_order(root)]
    n = rand(available)
    target::T = find_num(root, n)
    while target.nchild == 0 || any([ch.nchild == 0 for ch in target.children])
        n = rand(available)
        target = find_num(root, n)
    end

    # proportion of slide move is randomly selected
    proportion::Float64 = rand()
    # pick a random child
    child::T = rand(target.children)

    # calculate and set new values
    move!(target, child, proportion)

end # function slide!

"""
slide(root::T)::T where T<:GeneralNode

This functin performs a slide move on an intermediate node. The node is moved
upwards or downwards on the path specified by its mother and one of its
daughters.

Returns root Node of new tree.

* `root` : root Node of tree.
"""
function slide(root::T)::T where {T<:GeneralNode}
    new_root = deepcopy(root)
    slide!(new_root)
    return new_root
end


"""
swing!(root::T) where T<:GeneralNode

This function performs a swing node. A random non-leave node is selected and
moved along the path specified by its two children.

* `root` : root Node of tree.
"""
function swing!(root::T) where {T<:GeneralNode}

    available = [node.num for node in post_order(root)]
    n = rand(available)
    target::T = find_num(root, n)
    
    while target.nchild < 2
        n = rand(available)
        target = find_num(root, n)
    end
    proportion::Float64 = rand()

    child1 = target.children[1]
    child2 = target.children[2]

    # calculate and set new values
    move!(child1, child2, proportion)
end # function swing!

"""
swing(root::T)::T where T<:GeneralNode

This function performs a swing node. A random non-leave node is selected and
moved along the path specified by its two children.

Returns root Node of new tree.

* `root` : root Node of tree.
"""
function swing(root::T)::T where {T<:GeneralNode}
    new_root = deepcopy(root)
    swing!(new_root)
    return new_root
end

"""
    move!(node1::T, node2::T, proportion::Float64) where T <:GeneralNode

Change the incoming length of node1 and node2 while keeping their combined length
constant.

* `node1` : Node whose inc_length will be modified; this node's inc_length will be the total inc_length of both nodes, times proportion.

* `node2` : Node whose inc_length will be modified; this node's inc_length will be the remainder of total - the new inc_length value of node1.

* `proportion` : Float64, determines proportion of the inc_length of both nodes assigned to node1.
"""
function move!(node1::T, node2::T, proportion::Float64) where {T<:GeneralNode}
    total::Float64 = node1.inc_length + node2.inc_length
    fp::Float64 = total * proportion
    sp::Float64 = total - fp
    node1.inc_length = fp
    node2.inc_length = sp
end # function move!


### Altering the length of the tree ###

"""
    change_edge_length!(root::T) where T <:GeneralNode

Pick a random node and increase or decrease its length randomly.

* `root` : root node of tree.
"""
function change_edge_length!(root::T) where {T<:GeneralNode}
    available = [node.num for node in post_order(root)]
    n = rand(available)
    target::T = find_num(root, n)
    while target.root
        n = rand(available)
        target = find_num(root, n)
    end
    factor = abs(randn())
    target.inc_length *= factor
end


