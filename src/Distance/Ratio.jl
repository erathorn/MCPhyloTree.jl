#################### Ratio ####################
"""
    Ratio
    
This struct represents a ratio of vectors of tree edges. 

* `e_length` : Specifies the square root of the sum of the lengths of all e_edges.
* `f_length` : Specifies the square root of the sum of the lengths of all f_edges.
* `e_edges` : The edges of the first tree used in the geodesic calculation.
* `f_edges` : The edges of the second tree. 
"""
mutable struct Ratio{T<:GeneralNode}
    e_length::Float64
	f_length::Float64
	e_edges::Vector{T}
	f_edges::Vector{T}

    Ratio() = new{GeneralNode}(0.0, 0.0, [], [])
    Ratio(e::Vector{T}, f::Vector{T}) where T<:GeneralNode = new{typeof(e[1])}(geo_avg(e), geo_avg(f), e, f)
end # Ratio


"""
    add_e_edge!(ratio::Ratio, node::T) where T<:GeneralNode

--- INTERNAL ---
This function adds a new edge to one side of a ratio and adjusts the ratio length 
accordingly.

* `ratio` : ratio to which node is added
* `node` : node that is added to the ratio
"""
function add_e_edge!(ratio::Ratio, node::T) where T<:GeneralNode
    push!(ratio.e_edges, node)
    ratio.e_length = sqrt(ratio.e_length ^ 2 + node.inc_length ^ 2)
end # add_e_edge!


"""
    add_f_edge!(ratio::Ratio, node::T) where T<:GeneralNode

--- INTERNAL ---
See 'add_e_edge!'
"""
function add_f_edge!(ratio::Ratio, node::T) where T<:GeneralNode
    push!(ratio.f_edges, node)
    ratio.f_length = sqrt(ratio.f_length ^ 2 + node.inc_length ^ 2)
end # add_f_edge!


"""
    get_ratio(r::Ratio)::Float64

--- INTERNAL ---
Computes the ratio of a ratio struct.
"""
function get_ratio(r::Ratio)::Float64
    r.e_length / r.f_length
end # get_ratio


"""
combine(r1::Ratio, r2::Ratio)::Ratio

-- INTERNAL ---
Combine two ratios to make a new one.
"""
function combine(r1::Ratio, r2::Ratio)::Ratio
    r::Ratio = Ratio()
    edges::Vector{T where T<:GeneralNode} = []
    if length(r1.e_edges) == 0 && length(r2.e_edges) == 0
        r.e_length(sqrt(r1.e_length ^ 2 + r2.e_length ^ 2))
    else
        append!(edges, deepcopy(r1.e_edges), deepcopy(r2.e_edges))
        r.e_edges = edges
        r.e_length = geo_avg(r.e_edges) 
    end # if/else

    if length(r1.f_edges) == 0 && length(r2.f_edges) == 0
        r.f_length(sqrt(r1.f_length ^ 2 + r2.f_length ^ 2))
    else 
        append!(edges, deepcopy(r1.f_edges), deepcopy(r2.f_edges))
        r.f_edges = edges
        r.f_length = geo_avg(r.f_edges) 
    end # if/else
    return r
end # combine


#################### Ratio Sequence ####################
"""
    RatioSequence
    
This data type holds a sequence of ratios and a combine code.

* `ratios` : A vector of ratios, making up the sequence.
* `combine_code` : A code, relevant when RatioSequence needs to be combined.
"""
mutable struct RatioSequence
    ratios::Vector{Ratio}
    combine_code::Int64

    RatioSequence() = new([], 0)
end # RatioSequence

const EdgeLengths = Tuple{Float64, Float64}

"""
    get_distance(rs::RatioSequence)::Float64

--- INTERNAL ---
Returns the distance of a RatioSequence.
"""
function get_distance(rs::RatioSequence)::Float64
    distance²::Float64 = 0
    r::Ratio = Ratio()
    for i in 1:length(rs.ratios)
        r = rs.ratios[i]
        distance² += (r.e_length + r.f_length) ^ 2
    end # for
    return sqrt(distance²)
end # get_distance


"""
    interleave(rs1::RatioSequence, rs2::RatioSequence)::RatioSequence

--- INTERNAL ---
Interleaves the ratio sequences rs1 and rs2 after combining them to get the ascending ratio
sequence with the min distance. Returns a new ratio sequence.
"""
function interleave(rs1::RatioSequence, rs2::RatioSequence)::RatioSequence
    combined1::RatioSequence = get_non_desc_rs_with_min_dist(rs1)
    combined2::RatioSequence = get_non_desc_rs_with_min_dist(rs2)

    interleaved_rs = RatioSequence()
    ind1::Int64 = 1
    ind2::Int64 = 1
    while (ind1 <= length(combined1.ratios) && ind2 <= length(combined2.ratios))
        if get_ratio(combined1.ratios[ind1]) <= get_ratio(combined2.ratios[ind2])
            push!(interleaved_rs.ratios, combined1.ratios[ind1])
            ind1 += 1
        else
            push!(interleaved_rs.ratios, combined2.ratios[ind2])
            ind2 += 1
        end # if/else
    end # while

    while ind1 <= length(combined1.ratios)
        push!(interleaved_rs.ratios, combined1.ratios[ind1])
        ind1 += 1
    end # while

    while ind2 <= length(combined2.ratios)
        push!(interleaved_rs.ratios, combined2.ratios[ind2])
        ind2 += 1
    end # while

    return interleaved_rs
end # interleave


"""
    get_non_desc_rs_with_min_dist(rs::RatioSequence)::RatioSequence

--- INTERNAL ---
This function returns a combined RatioSequence, so that the ratios are non-descending.
"""
function get_non_desc_rs_with_min_dist(rs::RatioSequence)::RatioSequence
    length(rs.ratios) < 2 && return rs
    combined_rs::RatioSequence = deepcopy(rs)
    i::Int64 = 1
    combine_code::Int64 = 0
    combined_ratio::Ratio = Ratio()
    cc_array::Vector{Int64} = zeros(Int64, length(rs.ratios) - 1) .+ 2
    a::Int64 = 1

    while i < length(combined_rs.ratios) - 1
        if get_ratio(combined_rs.ratios[i]) > get_ratio(combined_rs.ratios[i+1])
            combined_ratio = combine(combined_rs.ratios[i], combined_rs.ratios[i+1])
            deleteat!(combined_rs.ratios,[i, i+1])
            insert!(combined_rs.ratios, i, combined_ratio)
            cc_array[a] = 1
            if i > 1
                i -= 1
                while cc_array[a] == 1
                    a -= 1
                end # while
            else
                while (a < length(rs.ratios) - 1) && cc_array[a] != 2
                    a += 1
                end # while
            end # if/else
        else
            cc_array[a] = 0
            i += 1
            while (a < length(rs.ratios) - 1) && cc_array[a] != 2
                a += 1
            end # while
        end # if/else
    end # while
    
    for k in 1:(length(rs.ratios) - 1)
        if cc_array[k] == 1
            combine_code = combine_code + 2 ^ k
        end # if
    end # for
    combined_rs.combine_code = combine_code
    return combined_rs
end # get_non_desc_rs_with_min_dist