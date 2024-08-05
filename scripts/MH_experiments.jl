import Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
Pkg.instantiate()

using BalancedMinimumEvolution
import BalancedMinimumEvolution as BME
using CSV, DataFrames, Statistics, JuMP, Combinatorics
import Random

# Matheuristic vs neighborhood joining 
datasets = ["01-Primates12", "02-M17", "03-M18", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "RDSM32", "RDSM64"]
begin
    experiments = [(id = "c"*string(Int(e))*"r"*string(Int(r)) , complete = e, nj_criterion = r) for (e,r) in Iterators.product((true, false),(true, false))]
    if !isfile("data/MH_results.csv") 
        df = DataFrame(dataset=[], fastmemh = [], fastme_b = [], fastme_i = [], fastme_o = [], fastme_n = [], fastme_u = [], heuristic = [], gap = [], nj = [], complete = [], nj_criterion = [], time = [], radius_complete = [], radius_relaxation= [],lp_gap=[])
        CSV.write("data/MH_results.csv", df)
    end
    global counter = 0
    df = CSV.read("data/MH_results.csv", DataFrame)
    for dataset in datasets
        for experiment in experiments
            global counter += 1
            if nrow(df) >= counter 
                continue #checkpointing
            end
            experiment_id = experiment.id
            complete = experiment.complete
            nj_criterion = experiment.nj_criterion
            println("-"^90)
            println(dataset)
            println(experiment_id)
            path = joinpath("data", dataset, dataset*".txt")
            exp_path = joinpath("data",dataset,"experiment$experiment_id")
            if !isdir(exp_path)
                mkdir(joinpath("data",dataset,"experiment$experiment_id"))
            end
            tree_path = joinpath(exp_path, "tree.nwk")
            D = read_distance_matrix(path);
            n = only(unique(size(D)))
            # Repair
            g = Graph(n+1)
            c = n+1
            for i in 1:c-1
                add_edge!(g, i, c)
            end
            D_ = BME.extend_distance(D);

            rtime = @elapsed gmh = BME.LP_heuristic(g, D_, c, complete = complete, nj_criterion = nj_criterion)
            tl = BME.tree_length(path_length_matrix(gmh), D)

            gfmmh = try
               open("tmp_tree.nwk", "w") do f
                    println(f, ubt_to_nwk(gmh))
                    fastme_local_search(path, "tmp_tree.nwk", inittree = true)
                end
            finally
                rm("tmp_tree.nwk")            
            end
            tlfmmh = tree_length(path_length_matrix(gfmmh), D)

            tlfastme = Dict{String, Float64}()
            for method in ["b", "i", "o", "n", "u"]
                try
                    gfm = fastme_local_search(path, tree_path, methods = [method])
                catch e
                    tlfastme[method] = Inf
                end

                tlfastme[method] = tree_length(path_length_matrix(gfm), D)
            end

            gnj = neighbor_joining(D)
            tlnj = tree_length(path_length_matrix(gnj), D)

            model = BME.BMEP_MILP(g, D_, c; relax = true)
            optimize!(model)
            τ_tilde = similar(D)
            for ((i,j), dist) in value.(model[:τ]).data 
                τ_tilde[i,j] = dist
                τ_tilde[i,i] = 0.
                τ_tilde[j,j] = 0.
            end
            radius_complete = maximum(abs.(path_length_matrix(gmh) .- τ_tilde[1:n,1:n]))
            
            model = BME.BMEP_MILP(g, D_, c; relax = true, complete = false)
            optimize!(model)
            for ((i,j), dist) in value.(model[:τ]).data 
                τ_tilde[i,j] = dist
                τ_tilde[i,i] = 0.
                τ_tilde[j,j] = 0.
            end
            radius_relaxation = maximum(abs.(path_length_matrix(gmh) .- τ_tilde[1:n,1:n]))
            lp_gap = tl/tree_length(τ_tilde, D)-1

            push!(df,(dataset=dataset, fastmemh = tlfmmh, fastme_b = tlfastme["b"], fastme_i = tlfastme["i"], fastme_o = tlfastme["o"], fastme_n = tlfastme["n"], fastme_u = tlfastme["u"], heuristic = tl, gap = tl/tlfastme-1, nj = tlnj, complete = complete, nj_criterion = nj_criterion, time = rtime, radius_complete=radius_complete, radius_relaxation=radius_relaxation, lp_gap = lp_gap), promote = true)
            CSV.write("data/MH_results.csv",df)            
        end
    end
end
