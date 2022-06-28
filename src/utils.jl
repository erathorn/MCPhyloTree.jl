"""
    lcp(str1::T, str2::T)::T where T <: AbstractString

Get the longest common prefix.
"""
function lcp(str1::T, str2::T)::T where T <: AbstractString
    minl = min(length(str1), length(str2))
    outs::T = ""
    minl == 0 && return outs
    @inbounds for i in 1:minl
      str1[i] == str2[i] ? outs *= str1[i] : return outs
    end # for
    return outs
end
