using BalancedMinimumEvolution
import BalancedMinimumEvolution as BME
using Graphs, JuMP

datasets = ["01-Primates12", "02-M17", "03-M18", "04-SeedPlants500", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "09-M82"]
dataset = datasets[1]
path = joinpath("data", dataset, dataset*".txt")
tree_path = path*"_tree.nwk"
D = read_distance_matrix(path);
n = only(unique(size(D)))

gfm = fastme_local_search(path, tree_path, inittree = false)
plmfm = path_length_matrix(gfm)
@show tlfastme = tree_length(plmfm, D)

g, c = star_graph(D)

D_ = BME.extend_distance(D);
τ = similar(D_);

model = BME.BMEP_MILP(g, D_ ,c, relax = false)
set_attribute(model, "TimeLimit", 3600)
# for i in 1:n
#     for j in (i+1):n
#         set_start_value(model[:_x][i,j,plmfm[i,j]], 1)
#     end
# end

# set_attribute(model, MOI.NumberOfThreads(), Threads.nthreads())
optimize!(model) 

for ((i,j), dist) in value.(model[:τ]).data 
    τ[i,j] = dist
    τ[i,i] = 0.
end

τ = τ[1:n,1:n]
τfm = path_length_matrix(gfm)
τdiff = τ .- τfm

tlmip = tree_length(τ, D)
tlmip/tlfastme

gopt = UBT_from_PLM(τ)

@assert τ == path_length_matrix(gopt)
path_length_matrix(gopt) .- τ

include("../src/plotting.jl")
plot_ubt(gopt)

maximum(τ)
maximum(τfm)
extrema(τdiff)
count(!=(0), τdiff)/2 / (n*(n-1))/2
