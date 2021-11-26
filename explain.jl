using Revise
using Pkg
Pkg.activate(".")
using MCPhyloTree

# trees = ParseNewick("Drav_mytrees_1.nwk")
trees = ParseNewick("gs_tree_JC_20-60-letters.nwk")


t1 = deepcopy(trees[1])
t2 = deepcopy(trees[2])
MCPhyloTree.number_nodes!(t1)
MCPhyloTree.number_nodes!(t2)
MCPhyloTree.set_binary!(t1)
MCPhyloTree.set_binary!(t2)

non_common_nodes = MCPhyloTree.splitOnCommonEdge(t1, t2)