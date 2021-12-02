using Revise
using Pkg
Pkg.activate(".")
using MCPhyloTree

# trees = ParseNewick("Drav_mytrees_1.nwk")
trees = ParseNewick("gs_tree_JC_20-60-letters.nwk")

t1 = deepcopy(trees[1])
t2 = deepcopy(trees[2])
MCPhyloTree.initialize_tree!(t1; height=false)
MCPhyloTree.initialize_tree!(t2; height=false)

cn = MCPhyloTree.getCommonEdges(t1, t2)

non_common_nodes = MCPhyloTree.splitOnCommonEdge(t1, t2)




