"""
    NNI(root::T, target::T, lor::Bool)::Int64   where T<:AbstractNode

This function does a nearest neighbour interchange (NNI) move on the tree specified
by `root`. The parameter `target` specifies the node which performs the interchange
move using the left or right child of the target node. If the left child should
be used `lor=true`.

The function returns 1 if the move was successful and 0 else.

* `root` : root node of tree on which to perform the NNI.

* `target` : specific node of tree to interchange.

* `lor` : Bool; "true" uses the left child of `target,` "false," the right child.
"""
function NNI!(root::T, target::T, lor::Bool)::Int64  where T<:AbstractNode
    # NNI move would be illegal
    if target.nchild == 0 || isroot(target)
        return 0
    end # if

    ancestor::T = get_mother(target)
    sister::T = get_sister(target)

    ychild::T = remove_child!(target, lor)
    xchild::T = remove_child!(ancestor, sister)

    add_child!(target, sister)
    add_child!(ancestor, ychild)

    set_binary!(root)

    return 1

end # function

"""
    NNI!(root::T, target::Int64)::Int64  where T<:AbstractNode

This function does a nearest neighbour interchange (NNI) move on the tree specified
by `root`. The target is identified by the number of the target node.

The function returns 1 if the move was successful and 0 else.

* `root` : root node of tree on which to perform the NNI.

* `target` : specific node of tree to interchange.
"""
function NNI!(root::T, target::Int64, lor::Bool)::Int64  where T<:AbstractNode
   tn::T = find_num(root, target)
   NNI!(root, tn, lor)
end #function


"""
    NNI!(root::T, target::Int64)::Int64  where T<:AbstractNode

This function does a nearest neighbour interchange (NNI) move on the tree specified
by `root`. The target is identified by the number of the target node.
The function returns 1 if the move was successfull and 0 else.
"""
function NNI!(root::T, target::Int64)::Int64  where T<:AbstractNode
    NNI!(Random.GLOBAL_RNG, root, target)
end

function NNI!(rng::Random.AbstractRNG, root::T, target::Int64)::Int64  where T<:AbstractNode
   tn::T = find_num(root, target)
   lor::Bool = 0.5 > rand(rng)
   NNI!(root, tn, lor)
end #function

"""
    NNI!(root::T)::Int64  where T<:AbstractNode

This function does a nearest neighbour interchange (NNI) move on the tree specified
by `root`. The target is identified by the number of the target node.

The function returns 1 if the move was successful and 0 else.

* `root` : root node of tree on which to perform the NNI.

"""
function NNI!(root::T)::Int64  where T<:AbstractNode
    NNI!(Random.GLOBAL_RNG, root)
end

function NNI!(rng::Random.AbstractRNG, root::T)::Int64  where T<:AbstractNode
    n = rand(rng, 1:size(root)[1])
    tn::T = find_num(root, n)
    lor::Bool = 0.5 > rand(rng)
    NNI!(root, tn, lor)
end #function

"""
    NNI(root::T)::T  where T<:AbstractNode

This function does a nearest neighbour interchange (NNI) move on the tree specified
by `root.`

Returns a mutated copy while leaving the original tree intact.

* `root` : root node of tree on which to perform the NNI.

"""
function NNI(root::T)::T where T<:AbstractNode
    new_root = deepcopy(root)
    NNI!(new_root)
    return new_root
end


function NNI(lm_r::A) where A<:AbstractArray{<:Real, 2}
    lm = deepcopy(lm_r)
    NNI!(lm)
    return lm
end

function NNI!(lm_r::A) where A<:AbstractArray{<:Real, 2}
    NNI!(Random.GLOBAL_RNG,lm_r)
end

function NNI!(rng::Random.AbstractRNG, lm_r::A) where A<:AbstractArray{<:Real, 2}
    NNI!(lm_r, rand(rng, axes(lm_r, 2)))
end

function NNI_int!(lm_r, target, child)
    lm = hcat(lm_r, ones(size(lm_r, 1))) 
    _, mother = findmin(x->find_mother_helper(x, lm, target), axes(lm, 2))
    sister = findfirst(x -> lm[:, x] == (lm[:, mother] .!= lm[:, target]), axes(lm, 2))
    # remove child from target
    lm[:, target] .-= lm[:, child]
    # add sister to target
    lm[:, target] .+= lm[:, sister]
    
    lm_r .= lm[:, 1:end-1] 
    return 1
end

function NNI!(rng::Random.AbstractRNG, lm_r::A, target::Int) where A<:AbstractArray{<:Real, 2}
    if sum(lm_r[:, target]) == 1
        return 0
    end
    
    _, child1 = findmax(x->find_column_helper(x, lm_r, target), axes(lm_r, 2))
    child2 = findfirst(x -> lm_r[:, x] == (lm_r[:, child1] .!= lm_r[:, target]), axes(lm_r, 2))
    child = rand(rng) > 0.5 ? child1 : child2
    NNI_int!(lm_r, target, child)
end


function NNI!(lm_r::A, target::Int) where A<:AbstractArray{<:Real, 2}
    NNI!(Random.GLOBAL_RNG, lm_r, target)
end