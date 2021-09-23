# Moves

There are several operations to change the topology or any of the branch lengths of a tree. 
The different functions which are available in MCPhyloTree are listed here.

While the tree building methods in this package are generally built in a way that ensures 
that the nodes of the tree are fully initialized, i.e. they have a unique number, binary 
representation and a height, the move methods listed here do **not** automatically update
these fields. After using these methods you might want to run *initialize_tree!* or 
*update_tree!* to ensure the nodes remain / are initialized. 

**WARNING**: Do not run *initialize_tree!* after changing a tree, if you rely on the num 
field of the nodes to identify them. The same applies to the method *number_nodes* that is
called in *initalize_tree*.

## Edge Length

```@autodocs
Modules = [MCPhyloTree]
Pages = ["Moves/EdgeLength.jl"]
```

## NNI

```@autodocs
Modules = [MCPhyloTree]
Pages = ["Moves/NNI.jl"]
```

## Randomize

```@autodocs
Modules = [MCPhyloTree]
Pages = ["Moves/Randomize.jl"]
```

## Rerooting

```@autodocs
Modules = [MCPhyloTree]
Pages = ["Moves/Rerooting.jl"]
```

## SPR

```@autodocs
Modules = [MCPhyloTree]
Pages = ["Moves/SPR.jl"]
```
