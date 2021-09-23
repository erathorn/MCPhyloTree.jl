# Basics

This section provides an overview over the basic functionalities offered.

When building / modifying trees with methods like *add_child!* or *insert_node!*, it might
be necessarry to run *initialize_tree!* or *update_tree!* to ensure the nodes remain /are 
initialized. 

**WARNING**: Do not run *initialize_tree!* after changing a tree, if you rely on the num
field of the nodes to identify them. The same applies to the method *number_nodes* that is
called in *initalize_tree*.
## Basic Tree Functionalities

```@autodocs
Modules = [MCPhyloTree]
Pages   = ["Basics/Tree_Basics.jl"]
Filter =
```

## Tree Pruning

```@autodocs
Modules = [MCPhyloTree]
Pages = ["Basics/Tree_Pruning.jl"]
```

## Tree Search

```@autodocs
Modules = [MCPhyloTree]
Pages = ["Basics/Tree_Search.jl"]
```

## Tree Traversal

```@autodocs
Modules = [MCPhyloTree]
Pages = ["Basics/Tree_Traversal.jl"]
```

## Tree Representations

```@docs
newick(root::T) where T<:GeneralNode
```
