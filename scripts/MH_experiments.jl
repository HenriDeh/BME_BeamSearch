import Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
Pkg.instantiate()

using BalancedMinimumEvolution
import BalancedMinimumEvolution as BME
using CSV, DataFrames, Statistics, JuMP, Combinatorics
import Random

# Matheuristic vs neighborhood joining 
datasets = [["01-Primates12", "02-M17", "03-M18", "04-SeedPlants500", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "09-M82"];["100_rdpii_F81", "100_rdpii_F84", "100_rdpii_K2P", "100_rdpii_JC69"]]#,"200_rdpii_F81", "200_rdpii_F84", "200_rdpii_K2P", "200_rdpii_JC69","300_zilla_F81", "300_zilla_F84", "300_zilla_K2P", "300_zilla_JC69"]]; 
#datasets = [["RDSM32"]; ["RDSM64"];["RDSM128"];["RDSM256"];["RDSM512"];["RDSM1024"];["RDSM2048"]]
#datasets = ["05-M43"]
output_path = "data/MH_results.csv"
begin
    if !isfile(output_path)
        df = DataFrame(dataset= [], 
                        n = [],
                        method = [],
                        length = [],
                        time = [],
                        id = []
                        )
        CSV.write(output_path, df)
    end
    df = CSV.read(output_path, DataFrame)
    checkpoints = Set(df.id) 
    for dataset in datasets
        println("-"^90)
        println(dataset)
        path = joinpath("data", dataset, dataset*".txt")
        Dfull = read_distance_matrix(path);
        N = only(unique(size(Dfull)))
        maxn = min(45, N)
        minn = maxn
        while minn > 20 
            minn -= 10 
        end
        exp_path = joinpath("data",dataset,"matheuristic_outputs")
        if !isdir(exp_path)
            mkdir(exp_path)
        end
        for n in minn:10:maxn
            println("n = ",n)
            D = Dfull[1:n,1:n]
            for max_coi_cuts in [0,5,10]
                method = "LPNJ"*string(max_coi_cuts)
                id = dataset*"("*string(n)*")"*method
                if id in checkpoints
                    continue
                end
                print(method)
                tree_path = joinpath(exp_path, "tree.nwk")
                rtime = @elapsed gmh = BME.LP_heuristic(D, triangular = true, buneman = false, max_coi_cuts = max_coi_cuts, scale = false)
                tl = BME.tree_length(path_length_matrix(gmh), D)
                push!(df, (dataset = dataset, n = n, method = method, length = tl, time = rtime, id = id), promote = true)
                print(" ", tl, " SPR => ")
                open(tree_path, "w") do f
                    write(f, BME.ubt_to_nwk(gmh))
                end
                rtime += @elapsed gfm = fastme_local_search(D, tree_path, inittree = true)
                tl = BME.tree_length(path_length_matrix(gfm), D)
                println(tl)
                method = method*"-spr"
                id = dataset*"("*string(n)*")"*method
                push!(df, (dataset = dataset, n = n, method = method, length = tl, time = rtime, id = id), promote = true)
            end
            tlfastme = Dict{String, Float64}()
            methods = ["b", "i", "o", "n", "u"]
            for spr in (true, false)
                for method in methods
                    label = method*(spr ? "-spr" : "")
                    id = dataset*"("*string(n)*")"*label
                    if id in checkpoints
                        continue
                    end
                    fmtime = @elapsed begin
                        gfm = try
                            gfm = fastme_local_search(D, "tmp_tree.nwk", methods = [method], spr = spr)
                            tlfastme[label] = tree_length(path_length_matrix(gfm), D)
                            gfm
                        catch e
                            tlfastme[label] = Inf
                        finally
                            rm("tmp_tree.nwk", force = true)
                        end
                    end
                    push!(df, (dataset = dataset, n = n, method = label, length = tlfastme[label], time = fmtime, id = id), promote = true)
                end
            end
        end
        CSV.write(output_path,df)        
        println()    
    end
end
