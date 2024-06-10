using BalancedMinimumEvolution
using CSV, DataFrames, DataFramesMeta
include("../src/plotting.jl")
# using CairoMakie

df = CSV.read("data/LNS_results.csv", DataFrame)
show(df, allrows = true)
begin
    CSV.write("LNS.csv", DataFrame()) 
    tables_gaps = @groupby df :dataset 
    methods = ["iLPNJ-co","iLPNJ-re","iLPC-co","iLPC-re"]
    for tab in tables_gaps
        if all(==(0), tab.gap)
            continue
        end
        d = DataFrame(method=[], gap=[], dia = [], changes=[], ext = [])
        ds = first(tab.dataset)
        if startswith(ds, "0")
            ds = ds[4:end]
        end
        open("LNS.csv", "a") do f
            write(f, "\\midrule\n")
            write(f, "\\multicolumn{5}{c}{\\bf{$ds}}\\\\\n")
        end
        push!(d, ["FastME", round(first(tab.fastme),digits = 6), first(tab.diameter_fastme), "", ""])
        for i in 1:4
            push!(d, [methods[i] , string(abs(tab.gap[i])>=1e-5 ? round(tab.gap[i]*100,digits =3) : repr(round(tab.gap[i]*100, sigdigits=3)))*"\\%", tab.diameter_heuristic[i], string(Int(round(tab.proportion[i]*100,digits=0)))*"\\%", "\$"*join(split(tab.extrema[i][2:end-1], ","),"\\triangleright")*"\$"], promote = true)
        end
        CSV.write("LNS.csv", d, append = true, delim = '&', newline = "\\\\\n")
    end
end

table_gaps = @by df :dataset begin
    :FastME = round(first(:fastme), digits = 6)
    :Diameter = first(:diameter_fastme)
    $"-co-nj" = (:gap[i] == 0 ? "0\\%" : string(round(:gap[i]*100, digits = 3)) *"\\%")
    $"-co-nj_d" = :diameter_heuristic[i] 
    $"-co-nj_changes" = string(:pair_changed[i])*" ("*string(Int(round(:proportion[i]*100, digits = 0)))*"\\%)"
    $"-co-nj_ext" = join(split(:extrema[i][2:end-1], ","), " -")
    $"-re-nj" = (:gap[2] == 0 ? "0\\%" : string(round(:gap[2]*100, digits = 3)) *"\\%") 
    $"-re-nj_d" = :diameter_heuristic[2]
    $"-re-nj_changes" = string(:pair_changed[2])*" ("*string(Int(round(:proportion[2]*100, digits = 0)))*"\\%)"
    $"-re-nj_ext" = join(split(:extrema[2][2:end-1], ","), " -")
    $"-co" = (:gap[3] == 0 ? "0\\%" : string(round(:gap[3]*100, digits = 3)) *"\\%") 
    $"-co_d" = :diameter_heuristic[3]
    $"-co_changes" = string(:pair_changed[3])*" ("*string(Int(round(:proportion[3]*100, digits = 0)))*"\\%)"
    $"-co_ext" = join(split(:extrema[3][2:end-1], ","), " -")
    $"-re" = (:gap[4] == 0 ? "0\\%" : string(round(:gap[4]*100, digits = 3)) *"\\%")
    $"-re_d" = :diameter_heuristic[4]
    $"-re_changes" = string(:pair_changed[4])*" ("*string(Int(round(:proportion[4]*100, digits = 0)))*"\\%)"
    $"-re_ext" = join(split(:extrema[4][2:end-1], ","), " -")
end
@rtransform! table_gaps begin
    :dataset = startswith(:dataset, "0") ? :dataset[4:end] : :dataset
end;
@rename! table_gaps :Dataset = :dataset
@rsubset! table_gaps begin
    any(g->g≠"0\\%", [$"-co-nj",$"-re-nj",$"-co",$"-re"])
end
CSV.write("LNS.csv", table_gaps, delim = '&', newline = "\\\\\n")



df2 = CSV.read("data/repair_results.csv", DataFrame)
show(df2, allrows = true)
table_repair = @by df2 :dataset begin
    :FastME = round(first(:fastme), digits = 6)
    $"-co-nj" = (:gap[i] == 0 ? "0\\%" : string(round(:gap[i]*100, digits = 3)) *"\\%")
    $"-co-nj time" = string(Int(round(:time[i], digits = 0)))
    $"-re-nj" = (:gap[2] == 0 ? "0\\%" : string(round(:gap[2]*100, digits = 3)) *"\\%") 
    $"-re-nj time" = string(Int(round(:time[2], digits = 0)))
    $"-co" = (:gap[3] == 0 ? "0\\%" : string(round(:gap[3]*100, digits = 3)) *"\\%") 
    $"-co time" = string(Int(round(:time[3], digits = 0)))
    $"-re" = (:gap[4] == 0 ? "0\\%" : string(round(:gap[4]*100, digits = 3)) *"\\%")
    $"-re time" = string(Int(round(:time[4], digits = 0)))
    $"NJ" = string(round((first(:nj)/first(:fastme)-1)*100, digits=3)) * "\\%"
end;
@rtransform! table_repair begin
    :dataset = startswith(:dataset, "0") ? :dataset[4:end] : :dataset
end;
@rename! table_repair :Dataset = :dataset
CSV.write("repair.csv", table_repair, delim = '&', newline = "\\\\\n")

#Compare M18
gfm = ubtgraph_from_nwk("data/03-M18/experimentc0r0/tree.nwk")
gmh = ubtgraph_from_nwk("data/03-M18/experimentc0r0/tree.nwk_1")

ffm = plot_ubt(gfm, title = "FastME")
display(GLMakie.Screen(),ffm)
fmh = plot_ubt(gmh, title = "Matheuristic")
display(GLMakie.Screen(),fmh)
jldsave("M18.jld2"; ffm, fmh)

## Seedplants 10 hours
cd(joinpath(@__DIR__, ".."))
dataset = "04-SeedPlants500"
path = joinpath("data", dataset, dataset*".txt")
D = read_distance_matrix(path);
g10 = ubtgraph_from_nwk("data/04-SeedPlants500/extended/tree.nwk_10499")
plmg10 = path_length_matrix(g10)
tl = tree_length(plmg10, D)
gfm04 = ubtgraph_from_nwk("data/04-SeedPlants500/extended/tree.nwk")
plmfm = path_length_matrix(gfm04)
tlfm = tree_length(plmfm, D)
gap = tl/tlfm -1

plmdiff = plmfm .- plmg10
extrema(plmdiff)
maximum(plmg10)
(count(>(0), plmdiff)÷2)/ ((size(D)[1]*size(D)[1]-1)÷2)
