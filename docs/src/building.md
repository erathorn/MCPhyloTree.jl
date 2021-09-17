# Tree Building & Transformations

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
