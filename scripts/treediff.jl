using BalancedMinimumEvolution
using BenchmarkTools
using JuMP, Combinatorics, LinearAlgebra, Graphs
using DataStructures

datasets = ["01-Primates12", "02-M17", "03-M18", "04-SeedPlants500", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "09-M82"]
dataset = datasets[3]
path = joinpath("data", dataset, dataset*".txt")
tree_path = path*"_tree.nwk"
D = read_distance_matrix(path);
n = only(unique(size(D)))
LNS_matheuristic(path, 15, tree_path, k=1, max_time = 0.25*60, max_it = Inf, fastme_it = Inf, inittree = false, lifted_tri = false, weak_bun = false, contraction = true, adaptive = false, repair_iterate = 1, complete = false, relax = true)

gfm = fastme_local_search(path, tree_path, inittree = false)
gmh = ubtgraph_from_nwk(tree_path*"_1")
plmfm = path_length_matrix(gfm)
plmmh = path_length_matrix(gmh)

tree_length(plmfm, D)
tree_length(plmmh, D)

sum(plmfm)÷2
sum(plmmh)÷2
(sum(plmfm)÷2)/(sum(plmmh)÷2)
extrema(plmfm)
extrema(plmmh)

diffplm = plmfm .- plmmh
sum(diffplm)
extrema(diffplm)
sort(counter(diffplm[idx] for idx in eachindex(diffplm) if diffplm[idx] != 0).map, rev = true)
plmfm[argmax(diffplm)]
plmmh[argmax(diffplm)]

changes = sum(1 for d in diffplm if d != 0)÷2 #number of pairs with new tau 
changes/length(combinations(1:n,2)) #proportion of pairs with a new tau

sort(counter(plmfm[idx] for idx in CartesianIndices(plmfm) if diffplm[idx] != 0).map)
sort(counter(plmmh[idx] for idx in CartesianIndices(plmfm) if diffplm[idx] != 0).map)
m = counter(plmfm[idx] => plmmh[idx] for idx in CartesianIndices(plmfm) if diffplm[idx] != 0).map
ms=sort(m, rev = true, by = x->m[x])
