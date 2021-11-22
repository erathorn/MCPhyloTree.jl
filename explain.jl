using Revise
using Pkg
Pkg.activate(".")
using MCPhyloTree

trees = ParseNewick("./Drav_mytrees_1.nwk")

t1 = deepycopy(trees[1])
t2 = deepcopy(trees[3500])

common_nodes = MCPhyloTree.splitOnCommonEdge(t1, t2, get_leaves(t1))