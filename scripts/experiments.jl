import Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
Pkg.instantiate()

using BalancedMinimumEvolution
import BalancedMinimumEvolution as BME
using CSV, DataFrames, Statistics, JuMP, Combinatorics
import Random

#datasets = [["01-Primates12", "02-M17", "03-M18", "04-SeedPlants500", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "09-M82"];["100_rdpii_F81", "100_rdpii_84", "100_rdpii_K2P", "100_rdpii_JC69","200_rdpii_F81", "200_rdpii_84", "200_rdpii_K2P", "200_rdpii_JC69","300_zilla_F81", "300_zilla_84", "300_zilla_K2P", "300_zilla_JC69"]]; ["RDSM32_"].*string.(1:50); ["RDSM64_"].*string.(1:50);["RDSM128_"].*string.(1:50);["RDSM256_"].*string.(1:50);["RDSM512_"].*string.(1:50);["RDSM1024_"].*string.(1:50);["RDSM2048_"].*string.(1:50)]
datasets = [["01-Primates12", "02-M17", "03-M18", "04-SeedPlants500", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "09-M82"];["100_rdpii_F81", "100_rdpii_84", "100_rdpii_K2P", "100_rdpii_JC69","200_rdpii_F81", "200_rdpii_84", "200_rdpii_K2P", "200_rdpii_JC69","300_zilla_F81", "300_zilla_84", "300_zilla_K2P", "300_zilla_JC69"]; ["RDSM32"]; ["RDSM64"];["RDSM128"];["RDSM256"];["RDSM512"];["RDSM1024"];["RDSM2048"]]

# LNS Matheuristic
begin
    experiments = [(id = "c"*string(Int(e))*"r"*string(Int(r))*"_$i", exact = e, nj_criterion = r, trial=i) for (e,r,i) in Iterators.product((true, false),(true, false),1:10)]
    if !isfile("data/LNS_results.csv") 
        df = DataFrame(dataset=[], fastme = [], heuristic = [], gap = [], exact = [], nj = [], time = [], tree_path = [], diameter_fastme=[],diameter_heuristic=[],pair_changed =[], proportion = [], extrema=[],avg_change=[],trial=[])
        CSV.write("data/LNS_results.csv", df)
    end
    global counter = 0
    df = CSV.read("data/LNS_results.csv", DataFrame)
    for dataset in datasets
        for experiment in experiments
            Random.seed!(counter)
            global counter += 1
            if nrow(df) >= counter 
                continue #checkpointing
            end
            experiment_id = experiment.id
            exact = experiment.exact
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
            K = min(exact ? 20 : 30, n-n%5)
            time = n < 20 ? 3 : n > 75 ? 60 : 15
            new_tree_path = LNS_matheuristic(path, K, tree_path, max_time = time*60, max_it = Inf, fastme_it = Inf, inittree = false, nj_criterion = nj_criterion, repair_iterate = 1, exact = exact, relax = true, )
            ## Tree differences
            gfm = ubtgraph_from_nwk(tree_path)
            gmh = ubtgraph_from_nwk(new_tree_path)
            plmfm = path_length_matrix(gfm)
            plmmh = path_length_matrix(gmh)
            tlfm=tree_length(plmfm, D)
            tlmh=tree_length(plmmh, D)
            @show diafm=maximum(plmfm)
            @show diamh=maximum(plmmh)
            diffplm = plmfm .- plmmh
            @show changes = sum([1 for d in diffplm if d != 0], init = 0)÷2 #number of pairs with new tau 
            @show avg = mean([d for d in diffplm if d != 0]) #number of pairs with new tau 
            @show prop = changes/length(combinations(1:n,2)) #proportion of pairs with a new tau
            @show ext=extrema(diffplm)
            push!(df,(dataset=dataset, fastme = tlfm, heuristic = tlmh, gap = tlmh/tlfm-1, exact = exact, nj = nj_criterion, time = time, tree_path = new_tree_path, diameter_fastme=diafm,diameter_heuristic=diamh,pair_changed =changes, proportion = prop, extrema=ext,avg_change=avg,trial=experiment.trial), promote = true)
            CSV.write("data/LNS_results.csv", df)
        end
    end
end

# Matheuristic vs neighborhood joining 
datasets = ["01-Primates12", "02-M17", "03-M18", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "RDSM32", "RDSM64"]
begin
    experiments = [(id = "c"*string(Int(e))*"r"*string(Int(r)) , exact = e, nj_criterion = r) for (e,r) in Iterators.product((true, false),(true, false))]
    if !isfile("data/repair_results.csv") 
        df = DataFrame(dataset=[], fastme = [], heuristic = [], gap = [], nj = [], exact = [], nj_criterion = [], time = [], radius_exact = [], radius_relaxation= [],lp_gap=[])
        CSV.write("data/repair_results.csv", df)
    end
    global counter = 0
    df = CSV.read("data/repair_results.csv", DataFrame)
    for dataset in datasets
        for experiment in experiments
            global counter += 1
            if nrow(df) >= counter 
                continue #checkpointing
            end
            experiment_id = experiment.id
            exact = experiment.exact
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

            rtime = @elapsed gmh = BME.LP_heuristic(g, D_, c, exact = exact, nj_criterion = nj_criterion)
            tl = BME.tree_length(path_length_matrix(gmh), D)

            gfm = ubtgraph_from_nwk(tree_path)
            tlfastme = tree_length(path_length_matrix(gfm), D)

            gnj = neighbor_joining(D)
            tlnj = tree_length(path_length_matrix(gnj), D)

            model = BME.MIP_complete(g, D_, c; relax = true)
            optimize!(model)
            τ_tilde = similar(D)
            for ((i,j), dist) in value.(model[:τ]).data 
                τ_tilde[i,j] = dist
                τ_tilde[i,i] = 0.
                τ_tilde[j,j] = 0.
            end
            radius_exact = maximum(abs.(path_length_matrix(gmh) .- τ_tilde[1:n,1:n]))
            
            model = BME.MIP_reduced(g, D_, c; relax = true)
            optimize!(model)
            for ((i,j), dist) in value.(model[:τ]).data 
                τ_tilde[i,j] = dist
                τ_tilde[i,i] = 0.
                τ_tilde[j,j] = 0.
            end
            radius_relaxation = maximum(abs.(path_length_matrix(gmh) .- τ_tilde[1:n,1:n]))
            lp_gap = tl/tree_length(τ_tilde, D)-1

            push!(df,(dataset=dataset, fastme = tlfastme, heuristic = tl, gap = tl/tlfastme-1, nj = tlnj, exact = exact, nj_criterion = nj_criterion, time = rtime, radius_exact=radius_exact, radius_relaxation=radius_relaxation, lp_gap = lp_gap), promote = true)
            CSV.write("data/repair_results.csv",df)            
        end
    end
end
