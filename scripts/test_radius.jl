using BalancedMinimumEvolution
import BalancedMinimumEvolution as BME
using BenchmarkTools
using JuMP, Combinatorics, LinearAlgebra, Graphs

datasets = ["01-Primates12", "02-M17", "03-M18", "04-SeedPlants500", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "09-M82"]
dataset = datasets[3]
path = joinpath("data", dataset, dataset*".txt")
tree_path = path*"_tree.nwk"
D = read_distance_matrix(path);
n = only(unique(size(D)))
contraction = true
iterate = n
begin
    g, c = star_graph(D)

    D_ = BME.extend_distance(D);
    τ = similar(D_);

    while degree(g, c) > 3
        model = BME.MIP_complete(g, D_ ,c, relax = false)
        #model = BME.MIP_reduced(g, D_ ,c, relax = true)
        #set_silent(model)
        set_optimizer_attribute(model, "TimeLimit", 3600)
        optimize!(model) 
        for ((i,j), dist) in value.(model[:τ]).data 
            τ[i,j] = dist
        end
        for _ in 1:n
            degree(g, c) > 3 || break
            (i,j) = contraction ? BME.contraction_select(g,c,τ) : BME.nj_select(g,c,τ)
            v = BME.join_neighbors!(g, c, i, j)
            BME.update_distances!(D_, g, c, v, i, j, method = contraction ? :contraction : :nj)
            BME.update_distances!(τ, g, c, v, i, j, method = contraction ? :contraction : :nj)
        end
    end
end
tl = BME.tree_length(path_length_matrix(g), D)
g2 = Graph(n+1)
c = n+1
for i in 1:c-1
    add_edge!(g2, i, c)
end
model = BME.MIP_complete(g2, D_ ,c, relax = true)
optimize!(model) 
τ_tilde = similar(D)
for ((i,j), dist) in value.(model[:τ]).data 
    if i <= n && j <= n
        τ_tilde[i,j] = dist
    end
end
for i in 1:n
    τ_tilde[i,i] = 0
end
τ_tilde
radius = maximum(abs.(path_length_matrix(g) .- τ_tilde))



gfm = fastme_local_search(path, tree_path, inittree = false)
tlfastme = tree_length(path_length_matrix(gfm), D)

tl/tlfastme