using BalancedMinimumEvolution
import BalancedMinimumEvolution as BME
using Graphs, JuMP

datasets = ["01-Primates12", "02-M17", "03-M18", "04-SeedPlants500", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "09-M82"]
dataset = datasets[3]
path = joinpath("data", dataset, dataset*".txt")
tree_path = path*"_tree.nwk"
D = read_distance_matrix(path);
n = only(unique(size(D)))

gfm = fastme_local_search(path, tree_path, inittree = false)
@show tlfastme = tree_length(path_length_matrix(gfm), D)

n = size(D,1)
g = Graph(n+1)
c = n+1
for i in 1:c-1
    add_edge!(g, i, c)
end

D_ = BME.extend_distance(D);
τ = similar(D_);

model = BME.MIP_complete(g, D_ ,c, relax = false)
set_attribute(model, "TimeLimit", 7200)
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

gopt = make_graph(τ)

@assert τ == path_length_matrix(gopt)
path_length_matrix(gopt) .- τ

include("../src/plotting.jl")
plot_ubt(gopt)

maximum(τ)
maximum(τfm)
extrema(τdiff)
count(!=(0), τdiff)/2 / (n*(n-1))/2
