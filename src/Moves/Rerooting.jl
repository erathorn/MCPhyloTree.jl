function reroot(root::T, new_root::T)::T where T<:GeneralNode

    new_tree = deepcopy(root)
    root_node = find_by_num(new_tree, new_root.num)
    mother = get_mother(root_node)

    recursive_invert(mother, root_node)

    root_node.root = true
    new_tree.root = false


    set_binary!(root_node)

    return root_node
end

# TODO (maybe): this only works when nodes have distinct names
function reroot(root::T, new_root::String)::T where T<:GeneralNode

    new_tree = deepcopy(root)
    root_node = find_by_name(new_tree, new_root)

    mother = get_mother(root_node)

    recursive_invert(mother, root_node)

    root_node.root = true
    new_tree.root = false


    set_binary!(root_node)

    return root_node
end


function recursive_invert(old_mother::T, old_daughter::T)::T where T
    if old_mother.root == true
        # arrived at the root
        od = remove_child!(old_mother, old_daughter)
        add_child!(od, old_mother)
        return od
    end
        od1 = recursive_invert(get_mother(old_mother), old_mother)
        od = remove_child!(od1, old_daughter)
        add_child!(od, od1)
        return od
end

