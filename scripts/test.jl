using BalancedMinimumEvolution
import BalancedMinimumEvolution as BME
using JuMP, Combinatorics, LinearAlgebra, Graphs

datasets = ["01-Primates12", "02-M17", "03-M18", "04-SeedPlants500", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "09-M82", "RDSM32","RDSM64","RDSM128","RDSM256","RDSM512","RDSM1024","RDSM2048"]
dataset = datasets[5]
#dataset = "100_rdpii_F81"
path = joinpath("data", dataset, dataset*".txt")
tree_path = path*"_tree.nwk"
D = read_distance_matrix(path);
n = only(unique(size(D)))
#FastME
g = fastme_local_search(path, tree_path)
tree_length(path_length_matrix(g), D)
#LNS
final = LNS_matheuristic(path, 40, tree_path, max_time = 60*60, max_it = Inf, inittree = false, nj_criterion = false, complete = false)
g2 = fastme_local_search(path, final, inittree = true)
tree_length(path_length_matrix(g2), D)
