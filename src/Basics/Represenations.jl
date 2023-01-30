
"""
    newick(root::T)::String  where T<:AbstractNode
Creates a newick representation of the tree.

Returns a properly formatted newick String.

* `node` : root node of tree used to create the newick string.
"""
function newick(root::T)::String  where T<:AbstractNode
    # get the newickstring
    newickstring = newick(root, "")

    # Some polishing.
    newickstring = chop(newickstring)
    newickstring = string(newickstring, ";")
    return newickstring
end

"""
    newick(root::T, newickstring::AbstractString) where T<:AbstractNode

Do the newick recursion. This is meant as the internal iterator function.
"""
function newick(root::T, newickstring::AbstractString) where T<:AbstractNode
    if root.nchild != 0
        # internal node
        newickstring = string(newickstring, "(")
        for child in root.children
            newickstring = string(newick(child,newickstring))
        end # for
        newickstring = chop(newickstring)
        return string(newickstring,")", root.name, ":", round(root.inc_length, digits=6),",")

    else
        # leave
        return string(newickstring, root.name, ":", root.inc_length, ",")
    end # if
end


