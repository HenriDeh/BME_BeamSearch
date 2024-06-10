using BalancedMinimumEvolution
using Graphs
include("../src/tree_traversal.jl")

g = Graph(11)
add_edge!(g, 1,2)
add_edge!(g, 2,4)
add_edge!(g, 4,3)
add_edge!(g, 4,5)
add_edge!(g, 6,2)
add_edge!(g, 6,7)
add_edge!(g, 7,9)
add_edge!(g, 9,8)
add_edge!(g, 7,10)
add_edge!(g, 6,11)

true_allorder = [11,6,2,1,1,1,2,4,3,3,3,4,5,5,5,4,2,6,7,9,8,8,8,9,9,7,10,10,10,7,6,11,11]

allorder(g, 11) == true_allorder

sts = traversal_subtrees(g, 30, allorder(g,11))
