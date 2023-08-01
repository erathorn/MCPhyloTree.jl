# Tree Building & Transformations

The tree building methods listed here ensure that the nodes of the tree(s) they build are
fully initialized, i.e. they have a unique number, binary representation and a height.
Therefore there is no need to run *initialize_tree!* or *update_tree!* after running them.

## Matrix Representation

```@autodocs
Modules = [MCPhyloTree]
Pages = ["Building/Tree2Matrix.jl"]
```

## Newick Parsing

```@autodocs
Modules = [MCPhyloTree]
Pages = ["Building/ParseNewick.jl"]
```

## Build Trees from Matrices

```@docs
MCPhyloTree.from_df
MCPhyloTree.create_tree_from_leaves
MCPhyloTree.cov2tree
MCPhyloTree.from_leave_incidence_matrix
```

## Tree Estimation from Matrices

```@docs
MCPhyloTree.neighbor_joining
MCPhyloTree.upgma
```

## Consensus Tree computation

```@docs
MCPhyloTree.majority_consensus_tree
MCPhyloTree.loose_consensus_tree
MCPhyloTree.greedy_consensus_tree
```

## Tree Ladderizing

```@autodocs
Modules = [MCPhyloTree]
Pages = ["Building/Tree_Ladderizing.jl"]
```
