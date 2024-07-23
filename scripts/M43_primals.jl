using BalancedMinimumEvolution
using JuMP, Combinatorics, LinearAlgebra, Graphs

dataset = "05-M43"
path = joinpath("data", dataset, dataset*".txt")
#dirpath = mkdir("data/M43primals")
for n in 13:43
    D = read_distance_matrix(path)[1:n,1:n];
    open("tmpM43.txt", "w") do f
        println(f, n)
        for (i,r) in enumerate(eachrow(D))
            println(f, "$i "*join(r,' '))
        end
    end
    tree_path = "tmpM43_tree.nwk"
    g = fastme_local_search("tmpM43.txt", tree_path)
    PLM = path_length_matrix(g)
    open(joinpath(dirpath, "M43_$n.dat"), "w") do f
        println(f, "PLM:[")
        for r in eachrow(PLM)
            println(f, join(r,' '))
        end
        println(f,"]")
    end
    rm("tmpM43.txt")
    rm(tree_path)
    rm("tmpM43.txt_fastme_stat.txt")
end