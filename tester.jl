using Pkg
Pkg.activate(".")
using Revise
using MCPhyloTree
using AbstractTrees

tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
po_list = post_order(tree)
length(po_list)