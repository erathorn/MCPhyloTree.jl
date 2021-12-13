
"""
function print_ascii(root::N) where N<:AbstractNode

Prints the ASCII representation of the tree.
Adapted for Julia from the ETE toolkit implementation.

* `root`: Root node of tree to print.
"""
function print_ascii(root::N) where N<:AbstractNode
lines,_ = ascii(root,"-")
println("\n"*join(lines,"\n"))
end #function


function ascii(root::N,char1="-") where N<:AbstractNode
    node_name = root.name
    LEN = max(3,length(node_name))
    PAD = " "^LEN
    PA = " "^(LEN-1)
    if root.nchild != 0
        mids = []
        result = []
        for c in root.children
            if length(root.children) == 1
                char2 = "-"
            elseif c.num == root.children[1].num
                char2 = "/"
            elseif c.num == root.children[end].num
                char2 = "\\"
            else
                char2 = "-"
            end #ifelse
            clines,mid = ascii(c,char2)


            push!(mids,mid+length(result))
            for x in clines
                push!(result,x)
            end #for
            push!(result,"")
        end #for
        pop!(result)
        (lo, hi, final) = (mids[1],mids[end],length(result))
        prefixes = []
        if lo+1 > 0
            for x in 1:(lo+1)
                push!(prefixes,PAD)
            end #for
        end #if
        if hi-lo-1 > 0
            for  x in 1:(hi-lo-1)
                push!(prefixes,PA*"|")
            end #for
        end #if
        if final-hi >0
            for x in 1:(final-hi)
                push!(prefixes,PAD)
            end #for
        end #if
        mid = (Integer(round((lo+hi)/2,RoundDown)))+1

        prefixes[mid] = char1 * ("-"^(LEN-2)) * (prefixes[mid][end])
        if length(prefixes) == length(result)
            result = [p*l for (p,l) in zip(prefixes,result)]
        else
            nuresult = []
            for x in 1:min(length(prefixes),length(result))
                push!(nuresult,prefixes[x]*result[x])
            end #for
            result = nuresult
        end #ifelse
        stem = result[mid]
        #added to handle spacing between a parent and child, each with a long name
        optional = ""
        if !(length(stem) < length(node_name)+2) && isletter(stem[(length(node_name)+2)])
            optional = "-"
        end #if

        result[mid] = stem[1] * node_name * optional * stem[(length(node_name)+2):end]

        return result,mid-1
    else
        return [char1 * "-" * node_name], 0
    end #if
end #function
