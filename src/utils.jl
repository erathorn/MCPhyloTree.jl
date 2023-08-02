"""
    lcp(str1::T, str2::T)::T where T <: AbstractString

Get the longest common prefix.
"""
function lcp(str1::T, str2::T)::T where T <: AbstractString
    minl = min(length(str1), length(str2))
    out_ind = 0
    minl == 0 && return outs
    @inbounds for i in 1:minl
      str1[i] == str2[i] ? out_ind += 1 : break
    end # for
    return str1[1:out_ind]
end


function find_column_helper(x_r, lm, current_ind, sister_ind = -1) 
    if x_r == sister_ind
        return -Inf
    end
    x = lm[:, x_r]
    rv = sum(x)
    if rv >= sum(lm[:, current_ind]) || length(intersect(findall(isone.(x)), findall(isone.(lm[:, current_ind])))) == 0
        rv = -Inf
    end
    return rv
end


function find_mother_helper(x_r, lm, current_ind) 
    x = lm[:, x_r]
    rv = sum(x)
    if rv <= sum(lm[:, current_ind]) || length(intersect(findall(isone.(x)), findall(isone.(lm[:, current_ind])))) == 0
        rv = Inf
    end
    return rv
end
