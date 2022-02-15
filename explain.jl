using Revise
using Pkg
Pkg.activate(".")
using MCPhyloTree

# trees = ParseNewick("Drav_mytrees_1.nwk")
trees = ParseNewick("gs_tree_JC_20-60-letters.nwk")

# t1 = deepcopy(trees[1])
# t2 = deepcopy(trees[2])
t1 = ParseNewick("(((((A:1,B:1):0.88,(C:1,D:1):1)J:0.47,E:1):0.73,F:1):0.83,G:1);")
t2 = ParseNewick("(((((C:0.2,D:1):0.5,A:1):0.15,E:1):0.87,(B:1,G:1):0.42):0.7,F:1);")

geo = geodesic(t1, t2; verbose=true)
println("Geodesic distance of the two trees is: $geo")