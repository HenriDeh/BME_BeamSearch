import Pkg; Pkg.activate("."); Pkg.instantiate()
using BalancedMinimumEvolution

dataset = "04-SeedPlants500"
exp_dir = joinpath("data", dataset, "extended")
if isdir(exp_dir)
    rm(exp_dir, recursive = true)    
end
mkdir(exp_dir)
path = joinpath("data", dataset, dataset*".txt")
tree_path = joinpath(exp_dir, "tree.nwk")
D = read_distance_matrix(path);
n = only(unique(size(D)))
check_radius = false
hybrid_matheuristic(path, 40, tree_path, max_time = 10*60*60, max_it = Inf, fastme_it = Inf, inittree = false, nj_criterion = false, repair_iterate = 1, exact = false, relax = true, nj_repair = false, check_radius = check_radius)
