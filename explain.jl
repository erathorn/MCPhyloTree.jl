using Revise
using Pkg
Pkg.activate(".")
using MCPhyloTree

trees = ParseNewick("./Drav_mytrees_1.nwk")
t1 = trees[1]
t2 = trees[420]

MCPhyloTree.splitOnCommonEdge(t1, t2, get_leaves(t1))