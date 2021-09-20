
"""
    newick(root::T)::String  where T<:GeneralNode
Creates a newick representation of the tree.

Returns a properly formatted newick String.

* `node` : root node of tree used to create the newick string.
"""
function newick(root::T)::String  where T<:GeneralNode
    # get the newickstring
    newickstring = newick(root, "")

    # Some polishing.
    newickstring = chop(newickstring)
    newickstring = string(newickstring, ";")
    return newickstring
end

"""
    newick(root::T, newickstring::AbstractString) where T<:GeneralNode

Do the newick recursion. This is meant as the internal iterator function.
"""
function newick(root::T, newickstring::AbstractString) where T<:GeneralNode
    if root.nchild != 0
        # internal node
        newickstring = string(newickstring, "(")
        for child in root.children
            newickstring = string(newick(child,newickstring))
        end # for
        newickstring = chop(newickstring)
        return string(newickstring,")", root.name, ":", root.inc_length,",")

    else
        # leave
        return string(newickstring, root.name, ":", root.inc_length, ",")
    end # if
end


