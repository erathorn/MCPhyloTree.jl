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

start = parsing_newick_string("(((((4:1,5:1):0.88,(3a:1,3b:1):1):0.47,2:1):0.73,1:1):0.83,6:1);")
target= parsing_newick_string("(((((3a:0.2,3b:1):0.5,4:1):0.15,2:1):0.87,(5:1,6:1):0.42):0.7,1:1);")
non_common_nodes = MCPhyloTree.splitOnCommonEdge(start, target)