"""
    is_valid_newick_string(newick::String)
This function checks if the given string is valid: is the brackets number matches and if the string ends with ";"
"""

function is_valid_newick_string(newick::String) # TODO: possibly not necessary; could be done as part of recursion, MAYBE
    # Step one: does the stripped string ends with ';'
    if endswith(strip(newick), ";")
        # Step two: check for the equal amount of brackets
        bracket_level = 0
        for letter in newick # => char is a type in Julia
            if letter == '('
                bracket_level += 1
            elseif letter == ')'
                bracket_level -= 1
            end # elseif
        end # for
        if bracket_level != 0
            return false
        end # if
    else # same level as the endswith statement
        return false
    end # else
    return true
end


"""
    parse_name_length(newick::String)

This function parses two optional elements of the tree, name and length. In case, 
when neither of this is provided, empty string and nothing are returned.
"""

function parse_name_length(newick::S) where S<:AbstractString

    if occursin(':', newick)
        name, len = split(strip(newick), ':')
        if name == ""
            name = "no_name"
        end # if
        if len == ""
            len = 1.0
        end # if
        if occursin(';', len)
            len = chop(len)
        end # if
        return string(name), parse(Float64, len)
    end # main if

    if length(newick) < 1
        return "no_name", 1.0
    else
        return string(newick), 1.0
    end # if-else
end # function

"""
    parsing_newick_string(newick::String)

In this function main parsing process happens, it uses recursive method to parse newick 
formated string.
"""

function parsing_newick_string(newick::A)::GeneralNode{Float64, Int64} where A<:AbstractString
    newick = replace(newick, " " => "")

    if newick[end] == ';' # no need for semicolon
        newick = chop(newick)
    end # if

    if newick[1] != ')' && newick[1] != '(' && occursin(r"^[a-zA-Z]|[0-9]*+[:]?[0-9]*", newick)
        leaf_node = Node()
        name, len = parse_name_length(newick)
        leaf_node.name = name
        leaf_node.inc_length = len
        return leaf_node
        # base case; only triggered at end of recursion OR if a single node-tree is input

    else
        current_node = Node()
        index = findlast(')', newick)[1]
        childrenstring = SubString(newick, 2, index - 1) # ... so that we can remove the superfluous parentheses here
        child_list = Sibling_parse(String(childrenstring))

        for x in child_list # recursion happens here
            add_child!(current_node, parsing_newick_string(childrenstring[x[1]:x[2]]))
        end # for

        child_list = []
        info_of_current_node = split(newick, ")") # info of current node should always follow the last ")"
        if lastindex(info_of_current_node) == 1
            name, length = parse_name_length(newick)
        else
            test = string(info_of_current_node[lastindex(info_of_current_node)])
            name, len = parse_name_length(string(test))
            current_node.name = name
            current_node.inc_length = len
        end # else
        return current_node
    end # recursion part
    throw("You left recursion somehow.")
end # function


function Sibling_parse(childrenstring::A)::Vector{Tuple{Int,Int}} where A<:AbstractString # returns list of children of a node
    child_list = Tuple{Int,Int}[]
    bracket_depth = 0
    start = 1
    
    for (ind, x) in enumerate(childrenstring * ",") # splits string identified above into a list, where each element corresponds to a child of current_node
        if x == ',' && bracket_depth == 0
            #push!(child_list, childrenstring[start:ind-1])
            push!(child_list, (start, ind-1))
            start = ind+1
            continue
        end # if
        if x == '('
            bracket_depth += 1
        end # if
        if x == ')'
            bracket_depth -= 1
        end # if
        
    end # for
    return child_list
end # function




"""
    ParseNewick(s::String)::Union{GeneralNode, Array{GeneralNode, 1}}

This function takes a string - either a filename or a newick string - and reads the file /
string to return an array of trees (represented as Node objects). The file should solely 
consist of newick tree representations, separated by line. The function checks for proper 
newick formatting, and will return an error if the string / file is incorrectly formatted.

Newick string input: Returns the root of the tree represented by the newick string.
Filename input: Returns an Array of Nodes; each Node is the root of the tree represented by
a newick string in the file.

* `s` : newick string or name of file containing newick strings to parse.
"""
function ParseNewick(s::String)::Union{GeneralNode, Array{GeneralNode,1}}
    
    # check if the input string is a newick string & parse + return it if so
    if is_valid_newick_string(s) 
        tree = parsing_newick_string(s)
        initialize_tree!(tree)
        return tree
    end # if
    
    list_of_trees::Vector{String} = []
    try
        list_of_trees = open(s, "r") do file
            readlines(file)
        end # do
    catch
        throw(ArgumentError("Input must be valid newick string or file name. Check path to file & check validity of newick string by running is_valid_newick_string(string)"))
    end # try/catch
    
    list_of_newicks = GeneralNode[]
    for content in list_of_trees
        if content == ""
            continue
        end # if
        if !is_valid_newick_string(content)
            throw(ArgumentError("$content is not correctly formatted!"))
        end # if
        tree = parsing_newick_string(string(content))
        initialize_tree!(tree)
        push!(list_of_newicks, tree)
    end # for
    list_of_newicks
end # ParseNewick