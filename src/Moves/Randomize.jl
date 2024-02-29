

"""
    randomize!(root::T, num::Int64=100)::Nothing where T <:AbstractNode

This function randomizes the tree by performing a number of nearest
neighbour interchange (NNI) moves and a randomization of the branch lengths.
The number of NNI moves is specified in the parameter num.

* `root` : root node of tree to be edited.

* `num` : amount of NNI moves to perform.
"""
function randomize!(root::T, num::Int64=100)::Nothing where T <:AbstractNode
    randomize!(Random.GLOBAL_RNG, root, num)
end

function randomize!(rng::Random.AbstractRNG, root::T, num::Int64=100)::Nothing where T <:AbstractNode
    n_nodes = size(root)[1]
    i = 0
    while i < num
        n = rand(rng, 1:n_nodes)
        NNI!(rng, root, n)
        i += 1
    end
    blv = rand(rng, n_nodes)
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
    return randomize(Random.GLOBAL_RNG, root, num)::T where T <:AbstractNode
end
function randomize(rng::Random.AbstractRNG, root::T, num::Int64=100)::T where T <:AbstractNode
    R = deepcopy(root)
    randomize!(rng, R, num)
    R
end

