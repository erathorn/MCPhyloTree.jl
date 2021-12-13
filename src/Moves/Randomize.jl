

"""
    randomize!(root::T, num::Int64=100)::Nothing where T <:AbstractNode

This function randomizes the tree by performing a number of nearest
neighbour interchange (NNI) moves and a randomization of the branch lengths.
The number of NNI moves is specified in the parameter num.

* `root` : root node of tree to be edited.

* `num` : amount of NNI moves to perform.
"""
function randomize!(root::T, num::Int64=100)::Nothing where T <:AbstractNode
    n_nodes = size(root)[1]
    i = 0
    while i < num
        n = rand(1:n_nodes)
        NNI!(root, n)
        i += 1
    end
    blv = rand(n_nodes)
    set_branchlength_vector!(root, blv)
end

"""
    randomize(root::T, num::Int64=100)::T where T <:AbstractNode

This function returns a randomized copy of the tree by performing a number of nearest
neighbour interchange (NNI) moves and a randomization of the branch lengths.
The number of NNI moves is specified in the parameter num. 

* `root` : root node of tree to be edited.

* `num` : amount of NNI moves to perform.
"""
function randomize(root::T, num::Int64=100)::T where T <:AbstractNode
    R = deepcopy(root)
    randomize!(R, num)
    R
end