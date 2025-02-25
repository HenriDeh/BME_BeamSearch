using BalancedMinimumEvolution
import BalancedMinimumEvolution as BME
using CSV, DataFrames, Statistics, JuMP, Combinatorics, PrettyTables, DataFramesMeta

df = CSV.read("data/experimental_results.csv", DataFrame)
df.best = zeros(Int, length(df.length))
df.spr_method = ifelse.(occursin.("spr", df.method), 1, 0)

lower_bounds = CSV.read("data/lower_bounds.csv", DataFrame)
df = leftjoin(df, lower_bounds, on = [:dataset, :n], makeunique = true)

@transform! df begin
    :gap = (:length .- :lb) ./ :lb
end

TOTAL_TIME = sum(df.time)/3600/2 # Divide by two because SPR double the construction time

dsA = Set(["01-Primates12", "02-M17", "03-M18", "04-SeedPlants500", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "09-M82"])
dsB = Set(["200_rdpii_F81", "200_rdpii_F84", "200_rdpii_K2P", "200_rdpii_JC69", "300_zilla_F81", "300_zilla_F84", "300_zilla_K2P", "300_zilla_JC69"])
dsC = Set(["RDSM"*string(n)*i for n in 10:5:60 for i in 'a':'j'])
dsD = Set(["RIM"*string(n)*i for n in 10:5:60 for i in 'a':'j'])

function identify_best_methods(df)
    grouped = groupby(df, [:dataset, :n, :spr_method])
    for group in grouped
        best = minimum(group.length)
        for row in eachrow(group)
            if row.length == best
                row.best = 1
            end
        end
    end
end

identify_best_methods(df)


dfA = filter(row -> row.dataset in dsA, df)
dfB = filter(row -> row.dataset in dsB, df)
dfC = filter(row -> row.dataset in dsC, df)
dfD = filter(row -> row.dataset in dsD, df)

function dataset_summary(df)
    if any(occursin.(r"RIM|RDSM", df.dataset))
        n_data = length(unique(df.dataset))
    else
        n_data = length(unique(df.dataset)) * length(unique(df.n))
    end
    summary = @by df :method begin
        :sum_best = sum(:best)
        :per_best = sum(:best)/n_data 
        :spr_method = only(unique(:spr_method))
        :mean_gap = mean(:gap)
    end
    sort!(summary, :per_best, rev = true)
    groupby(summary, :spr_method)
end

summaryA = dataset_summary(dfA)
summaryB = dataset_summary(dfB)
summaryC = dataset_summary(dfC)
summaryD = dataset_summary(dfD)
let
    dftab = DataFrame(method = [m for m in unique(df.method) if !contains(m, "-spr")])
    for spr in (false, true)
        for (name, summary) in [("dsA", summaryA), ("dsB", summaryB), ("dsC", summaryC), ("dsD", summaryD)]
            subdf = select(summary[1+spr], Not(:spr_method))
            subdf.method = replace.(subdf.method, r"-spr" => "")
            subdf.per_best = round.(subdf.per_best .* 100, digits=1)
            subdf.mean_gap = round.(subdf.mean_gap .* 100, digits=2)
            rename!(subdf, :per_best => "per_best_$(name[3])", :mean_gap => "mean_gap_$(name[3])")
            dftab = outerjoin(dftab, select(subdf, Not(:sum_best)), on = :method, makeunique = true)
        end
    end
    dftab.B = [!isnothing(match(r"(?-i)beam_(\d+)", m)) ? parse(Int, string(match(r"(?-i)beam_(\d+)", m).captures[1])) : 0 for m in dftab.method]
    select!(dftab, :method, :B, Not([:method, :B]))
    for row in eachrow(dftab)
        row.method = replace(row.method, Dict(
            r"(?-i)n$" => "NJ",
            r"(?-i)u$" => "Unweighted NJ",
            r"(?-i)b$" => "TaxAdd BME",
            r"(?-i)i$" => "BIONJ",
            r"(?-i)o$" => "TaxAdd OLS",
            r"(?-i)beam_(\d+)_nj" => s -> "BeamNJ",
            r"(?-i)beam_(\d+).*?_false" => s -> "BeamLPNJ-Tri",
            r"(?-i)beam_(\d+).*?_true" => s -> "BeamLPNJ+Tri"
        )...)
    end
    dftab = @orderby dftab begin
        (!contains).(:method, "Beam")
        -1 .* :B
        :method
    end
    dftab.B = string.(dftab.B)
    replace!(dftab.B, "0" => "-")
    open("beam_tables.tex", "w") do io
        println(io, "\\begin{tabular}{llrrrrrrrr|rrrrrrrr}")
        println(io, "\\toprule")
        println(io, "&& \\multicolumn{8}{c}{\\textbf{Agglomerative}} & \\multicolumn{8}{c}{\\textbf{Agglomerative + SPR}} \\\\")
        println(io, "\\cmidrule(lr){3-10} \\cmidrule(lr){11-18}")
        println(io, "&& \\multicolumn{2}{c}{Set A} & \\multicolumn{2}{c}{Set B}& \\multicolumn{2}{c}{Set C}& \\multicolumn{2}{c}{Set D} & \\multicolumn{2}{c}{Set A} & \\multicolumn{2}{c}{Set B}& \\multicolumn{2}{c}{Set C}& \\multicolumn{2}{c}{Set D} \\\\")
        println(io, "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10} \\cmidrule(lr){11-12} \\cmidrule(lr){13-14} \\cmidrule(lr){15-16} \\cmidrule(lr){17-18}")
        println(io, "\\textbf{Method} & \\textbf{B} & \\textbf{Best} & \\textbf{Gap} & \\textbf{Best} & \\textbf{Gap} & \\textbf{Best} & \\textbf{Gap} & \\textbf{Best} & \\textbf{Gap} & \\textbf{Best} & \\textbf{Gap} & \\textbf{Best} & \\textbf{Gap} & \\textbf{Best} & \\textbf{Gap} & \\textbf{Best} & \\textbf{Gap} \\\\")
        println(io, "\\midrule")
        best_values = [isodd(i) ? maximum(dftab[!, col]) :  minimum(dftab[!, col]) for (i,col) in enumerate(names(dftab)[3:end])]
        for row in eachrow(dftab)
            row_values = [row[col] for col in names(dftab)]
            formatted_row = [row_values[1], row_values[2]]  # First two columns
            for (i, val) in enumerate(row_values[3:end])
                if val == best_values[i]
                    push!(formatted_row, "\\textbf{$val}")
                else
                    push!(formatted_row, string(val))
                end
            end
            println(io, join(formatted_row, " & "), " \\\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end
end

# Generate table for Data Set A with each value of n as columns
begin
    dfA_gap = filter(dfA) do r 
        !contains(r.method, "-spr")
    end
    dfA_gap = @by dfA_gap [:method, :n] begin
        :gap =  round(median(:gap)*100, sigdigits = 2)
    end
    dfA_gap = unstack(dfA_gap, :n, :gap)

    dfA_gap.B = [!isnothing(match(r"(?-i)beam_(\d+)", m)) ? parse(Int, string(match(r"(?-i)beam_(\d+)", m).captures[1])) : 0 for m in dfA_gap.method]
    select!(dfA_gap, :method, :B, Not([:method, :B]))
    for row in eachrow(dfA_gap)
        row.method = replace(row.method, Dict(
            r"(?-i)n$" => "NJ",
            r"(?-i)u$" => "Unweighted NJ",
            r"(?-i)b$" => "TaxAdd BME",
            r"(?-i)i$" => "BIONJ",
            r"(?-i)o$" => "TaxAdd OLS",
            r"(?-i)beam_(\d+)_nj" => s -> "BeamNJ",
            r"(?-i)beam_(\d+).*?_false" => s -> "BeamLPNJ-Tri",
            r"(?-i)beam_(\d+).*?_true" => s -> "BeamLPNJ+Tri"
        )...)
    end
    dfA_gap = @orderby dfA_gap begin
        (!contains).(:method, "Beam")
        -1 .* :B
        :method
    end
    dfA_gap.B = string.(dfA_gap.B)
    replace!(dfA_gap.B, "0" => "-")
    open("beam_tables_gap_A.tex", "w") do io
        println(io, "\\begin{tabular}{llrrrrrrrrrrr}")
        println(io, "\\toprule")
        println(io, "&& \\multicolumn{11}{c}{\\textbf{Instance size (n)}} \\\\")
        println(io, "\\cmidrule(lr){3-13}")
        println(io, "\\textbf{Method} & \\textbf{B} & 10 & 15 & 20 & 25 & 30 & 35 & 40 & 45 & 50 & 55 & 60 \\\\")
        println(io, "\\midrule")
        for row in eachrow(dfA_gap)
            row_values = [row[col] for col in names(dfA_gap)]
            println(io, join(row_values, " & "), " \\\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end
end
# Same with SPR
begin
    dfA_gap = filter(dfA) do r 
        contains(r.method, "-spr")
    end
    dfA_gap = @by dfA_gap [:method, :n] begin
        :gap =  round(median(:gap)*100, sigdigits = 2)
    end
    dfA_gap = unstack(dfA_gap, :n, :gap)

    dfA_gap.B = [!isnothing(match(r"(?-i)beam_(\d+)", m)) ? parse(Int, string(match(r"(?-i)beam_(\d+)", m).captures[1])) : 0 for m in dfA_gap.method]
    select!(dfA_gap, :method, :B, Not([:method, :B]))
    for row in eachrow(dfA_gap)
        row.method = replace(row.method, Dict(
            r"(?-i)n$" => "NJ",
            r"(?-i)u$" => "Unweighted NJ",
            r"(?-i)b$" => "TaxAdd BME",
            r"(?-i)i$" => "BIONJ",
            r"(?-i)o$" => "TaxAdd OLS",
            r"(?-i)beam_(\d+)_nj" => s -> "BeamNJ",
            r"(?-i)beam_(\d+).*?_false" => s -> "BeamLPNJ-Tri",
            r"(?-i)beam_(\d+).*?_true" => s -> "BeamLPNJ+Tri"
        )...)
    end
    dfA_gap = @orderby dfA_gap begin
        (!contains).(:method, "Beam")
        -1 .* :B
        :method
    end
    dfA_gap.B = string.(dfA_gap.B)
    replace!(dfA_gap.B, "0" => "-")
    open("beam_tables_gap_A.tex", "w") do io
        println(io, "\\begin{tabular}{llrrrrrrrrrrr}")
        println(io, "\\toprule")
        println(io, "&& \\multicolumn{11}{c}{\\textbf{Instance size (n)}} \\\\")
        println(io, "\\cmidrule(lr){3-13}")
        println(io, "\\textbf{Method} & \\textbf{B} & 10 & 15 & 20 & 25 & 30 & 35 & 40 & 45 & 50 & 55 & 60 \\\\")
        println(io, "\\midrule")
        for row in eachrow(dfA_gap)
            row_values = [row[col] for col in names(dfA_gap)]
            println(io, join(row_values, " & "), " \\\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end
end

## Compute time for each method
begin
    dftime = filter(df) do r 
        !contains(r.method, "-spr")
    end

    dftime  = @by dftime [:method, :n] begin
        :time = ifelse(median(:time) < 0.01, round(median(:time), sigdigits = 2), round(median(:time), digits = 2))
    end
    dftime = unstack(dftime, :n, :time)

    dftime.B = [!isnothing(match(r"(?-i)beam_(\d+)", m)) ? parse(Int, string(match(r"(?-i)beam_(\d+)", m).captures[1])) : 0 for m in dftime.method]
    select!(dftime, :method, :B, Not([:method, :B]))
    for row in eachrow(dftime)
        row.method = replace(row.method, Dict(
            r"(?-i)n$" => "NJ",
            r"(?-i)u$" => "Unweighted NJ",
            r"(?-i)b$" => "TaxAdd BME",
            r"(?-i)i$" => "BIONJ",
            r"(?-i)o$" => "TaxAdd OLS",
            r"(?-i)beam_(\d+)_nj" => s -> "BeamNJ",
            r"(?-i)beam_(\d+).*?_false" => s -> "BeamLPNJ-Tri",
            r"(?-i)beam_(\d+).*?_true" => s -> "BeamLPNJ+Tri"
        )...)
    end
    dftime = @orderby dftime begin
        (!contains).(:method, "Beam")
        -1 .* :B
        :method
    end
    dftime.B = string.(dftime.B)
    replace!(dftime.B, "0" => "-")
    open("beam_tables_time.tex", "w") do io
        println(io, "\\begin{tabular}{llrrrrrrrrrrr}")
        println(io, "\\toprule")
        println(io, "&& \\multicolumn{11}{c}{\\textbf{Instance size (n)}} \\\\")
        println(io, "\\cmidrule(lr){3-13}")
        println(io, "\\textbf{Method} & \\textbf{B} & 10 & 15 & 20 & 25 & 30 & 35 & 40 & 45 & 50 & 55 & 60 \\\\")
        println(io, "\\midrule")
        for row in eachrow(dftime)
            row_values = [row[col] for col in names(dftime)]
            println(io, join(row_values, " & "), " \\\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end
end

