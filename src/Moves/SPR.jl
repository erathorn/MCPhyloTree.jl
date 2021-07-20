"""
    SPR(original_root::GeneralNode)::GeneralNode
Performs SPR on tree. Takes a copy of root of the tree;
Returns a copy of root of altered tree. Throws error if tree is improperly formatted.
"""
function SPR(original_root::GeneralNode)
    root = deepcopy(original_root)
    SPR!(root)
    return root
end #func

"""
        SPR!(root::GeneralNode)::GeneralNode
    Performs SPR on tree in place. Takes reference to root of tree;
    Returns reference to root of altered tree. Throws error if tree is improperly formatted.
"""
function SPR!(root::GeneralNode)
    if length(post_order(root)) <= 2
        error("The tree is too small for SPR")
    end #if
    check_binary(root) || throw(ArgumentError("Not yet implemented for not binary trees"))
    perform_spr(root)
    return root
end #function

"""
        risky_SPR(root::GeneralNode)::GeneralNode
    Performs SPR on tree in place. Takes reference to root of tree
    Returns copy of root of altered tree. Does not check for correct formatting of tree.
"""
function risky_SPR(original_root::GeneralNode)
    root = deepcopy(original_root)
    return risky_SPR!(root)
end #function


"""
        risky_SPR!(root::GeneralNode)::GeneralNode
Performs SPR on tree in place.

Returns reference to root of altered tree. Does not check for correct formatting of tree.

* `root` : root node of tree.
"""
function risky_SPR!(root::GeneralNode)
    return perform_spr(root)
end #func

"""
    perform_spr(root::T) where T <: GeneralNode
performs SPR on binary tree.

Returns root of tree post-SPR.

* `root` : Node of tree on which to perform SPR.
"""
function perform_spr(root::T) where T <: GeneralNode
    # find node to move
    available = [n.num for n in post_order(root)]
    n = rand(available)
    tn::T = find_num(root, n) #this is the root of the subtree which will be moved
    while tn.root || get_mother(tn).root
        n = rand(available)
        tn = find_num(root, n) #this is the root of the subtree which will be moved
    end # while
    tn_mother = get_mother(tn)
    tn_sister = get_sister(tn)
    tn_gm = get_mother(tn_mother)
    remove_child!(tn_gm, tn_mother)
    remove_child!(tn_mother, tn_sister)
    add_child!(tn_gm, tn_sister)
    # find target
    available = [n.num for n in post_order(root)]
    n = rand(available)
    target::T = find_num(root, n) #this is the target of the movement
    while target.root
        n = rand(available)
        target = find_num(root, n)
    end # while
    target_mother = get_mother(target)
    remove_child!(target_mother, target)
    add_child!(target_mother, tn_mother)
    add_child!(tn_mother, target)
    set_binary!(root)
    return root
end #func

"""
    perform_spr(root::T,subtree::T) where T <: GeneralNode
performs SPR on binary tree.

Returns root of tree post-SPR.

* `root` : Node of tree on which to perform SPR.
* `subtree`   : Root of subtree to be moved.
"""
function perform_spr(root::T,subtree::T) where T <: GeneralNode
    tn_mother = get_mother(subtree)
    tn_sister = get_sister(subtree)
    tn_gm = get_mother(tn_mother)
    remove_child!(tn_gm, tn_mother)
    remove_child!(tn_mother, tn_sister)
    add_child!(tn_gm, tn_sister)
    # find target
    available = [n.num for n in post_order(root)]
    n = rand(available)
    target::T = find_num(root, n) #this is the target of the movement
    while target.root
        n = rand(available)
        target = find_num(root, n)
    end # while
    target_mother = get_mother(target)
    remove_child!(target_mother, target)
    add_child!(target_mother, tn_mother)
    add_child!(tn_mother, target)
    set_binary!(root)
    return root
end #func

"""
    perform_spr(root::T,subtree::T,target::T) where T <: GeneralNode
performs SPR on binary tree.

Returns root of tree post-SPR.

* `root` : Node of tree on which to perform SPR.
* `subtree`   : Root of subtree to be moved.
* `target`    : Target of SPR.
"""
function perform_spr(root::T,subtree::T,target::T) where T <: GeneralNode
    @assert !subtree.root ["subtree cannot be the root node!"]
    @assert !get_mother(subtree).root ["subtree cannot be a child of the root node!"]
    tn_mother = get_mother(subtree)
    tn_sister = get_sister(subtree)
    tn_gm = get_mother(tn_mother)
    remove_child!(tn_gm, tn_mother)
    remove_child!(tn_mother, tn_sister)
    add_child!(tn_gm, tn_sister)
    @assert !target.root ["target cannot be the root!"]
    target_mother = get_mother(target)
    remove_child!(target_mother, target)
    add_child!(target_mother, tn_mother)
    add_child!(tn_mother, target)
    set_binary!(root)
    return root
end #func
