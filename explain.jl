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

geo = MCPhyloTree.geodesic(t1, t2)
dist = MCPhyloTree.get_distance(geo)
println("Geodesic distance of the two trees is: $dist")

