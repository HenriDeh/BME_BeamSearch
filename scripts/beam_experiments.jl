import Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
Pkg.instantiate()

using BalancedMinimumEvolution
import BalancedMinimumEvolution as BME
using CSV, DataFrames, Statistics, JuMP, Combinatorics, PrettyTables, DataFramesMeta
import Random

# beam search vs neighborhood joining 
datasets = ["01-Primates12", "02-M17", "03-M18", "04-SeedPlants500", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "09-M82"]
append!(datasets, [["200_rdpii_F81", "200_rdpii_F84", "200_rdpii_K2P", "200_rdpii_JC69"];["300_zilla_F81", "300_zilla_F84", "300_zilla_K2P", "300_zilla_JC69"]])
append!(datasets, ["RDSM"*string(n) for n in 10:5:60])
append!(datasets, ["RIM"*string(n) for n in 10:5:60])
output_path = "data/experimental_results.csv"
begin
    if !isfile(output_path)
        df = DataFrame(dataset = String[], n = Int[], method = String[], length = Float64[], time = Float64[], id = String[])
        CSV.write(output_path, df)
    end
    df = CSV.read(output_path, DataFrame)
    checkpoints = Set(df.id)
    for dataset in datasets
        is_biological_instance = !contains(dataset, r"(RIM|RDSM)")
        println("-"^90)
        println(dataset)
        path = joinpath("data", dataset, dataset*".txt")
        Dfull = read_distance_matrix(path);
        N = only(unique(size(Dfull)))
        if !is_biological_instance
            maxn = minn = N
        else
            maxn = 60
            minn = 10
        end
        exp_path = joinpath("data",dataset,"beam_outputs")
        if !isdir(exp_path)
            mkdir(exp_path)
        end
        for n in minn:10:maxn
            println("n = ",n)
            D = Dfull[1:n,1:n]
            for (tri, lp) in ((true, true), (false, true), (false, false))
                for B in [1,5,10,15]
                    method = "beam_"*string(B)*"_tri_"*string(tri)
                    id = dataset*"("*string(n)*")"*method
                    if id in checkpoints
                        continue
                    end
                    print(method)
                    tree_path = joinpath(exp_path, "tree.nwk")
                    timed_result = @timed trees = BME.beam_search(D, B, triangular = tri, buneman = false, scale = true, max_length = true, spr = false, lp = lp)
                    rtime = timed_result.time
                    tl = minimum(BME.tree_length(path_length_matrix(gmh), D) for gmh in trees)
                    push!(df, (dataset = dataset, n = n, method = method, length = tl, time = rtime, id = id), promote = true)
                    print(" ", tl, " SPR => ")
                    timed_result = @timed trees_spr = [fastme_local_search(D, g) for g in trees]
                    rtime += timed_result.time
                    tl = minimum(BME.tree_length(path_length_matrix(gfm), D) for gfm in trees_spr)
                    println(tl)
                    method = method*"-spr"
                    id = dataset*"("*string(n)*")"*method
                    push!(df, (dataset = dataset, n = n, method = method, length = tl, time = rtime, id = id), promote = true)
                end
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
                    timed_result = @timed begin
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
                    fmtime = timed_result.time
                    push!(df, (dataset = dataset, n = n, method = label, length = tlfastme[label], time = fmtime, id = id), promote = true)
                end
            end
        end
        CSV.write(output_path,df)        
        println()    
    end
end

begin
    open("beam_tables.md", "w") do f
        df = CSV.read(output_path, DataFrame)[:,1:4]
        gdf = groupby(df, :dataset)
        for sdf in gdf
            dataset = only(unique(sdf.dataset))
            println(f, "# $dataset")
            sdf1 = @rsubset sdf begin
                !contains(:method, r"-spr")
            end
            t1 = unstack(sdf1, :method, :n, :length)
            pretty_table(f, t1,backend = Val(:markdown), header = names(t1), highlighters = MarkdownHighlighter((data, i, j) -> j != 1 && (minimum(data[:, j]) == data[i,j]), MarkdownDecoration(bold = true)))
            println(f, "## + SPR")
            sdf2 = @rsubset sdf begin
                contains(:method, r"-spr")
            end
            t2 = unstack(sdf2, :method, :n, :length)
            pretty_table(f, unstack(sdf2, :method, :n, :length),backend = Val(:markdown), header = names(t2), highlighters = MarkdownHighlighter((data, i, j) -> j != 1 && (minimum(data[:, j]) == data[i,j]), MarkdownDecoration(bold = true)))
            println(f, "\n"*"-"^90)
        end
    end
end