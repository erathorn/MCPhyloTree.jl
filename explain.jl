using Revise
using Pkg
Pkg.activate(".")
using MCPhyloTree

trees = ParseNewick("Drav_mytrees_1.nwk")
# trees = ParseNewick("gs_tree_JC_20-60-letters.nwk")

t1 = deepcopy(trees[1])
t2 = deepcopy(trees[2])
MCPhyloTree.initialize_tree!(t1; height=false)
MCPhyloTree.initialize_tree!(t2; height=false)

cn = MCPhyloTree.get_common_edges(t1, t2)

nc = MCPhyloTree.split_on_common_edge(t1, t2)
leaves = get_leaves(nc[5][1])
leaf_dict = Dict(leaf.name => leaf.num for leaf in leaves)
internal_nodes1 = filter!(x -> x.nchild != 0, post_order(nc[5][1])[1:end-1])
internal_nodes2 = filter!(x -> x.nchild != 0, post_order(nc[5][2])[1:end-1])
mat = MCPhyloTree.get_incidence_matrix(internal_nodes1, internal_nodes2, leaf_dict)
bg = MCPhyloTree.build_bipartite_graph(mat, [n.inc_length for n in internal_nodes1], [n.inc_length for n in internal_nodes2])
MCPhyloTree.get_geodesic_nocommon_edges(t1, t2)