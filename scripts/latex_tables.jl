let
    df = CSV.read(output_path, DataFrame)[:,1:4]
    gdf = groupby(df, [:dataset, :n])
    cons = Dict(m => 0 for m in unique(df.method) if !occursin("-spr", m))
    n_beam_was_better = 0
    n_beam_was_worse = 0
    cons_spr = Dict(m => 0 for m in unique(df.method) if occursin("-spr", m))
    n_beam_spr_was_better = 0
    n_beam_spr_was_worse = 0
    n_data = length(gdf)
    println(n_data, " datasets")
    for sdf in gdf
        sdf1 = @rsubset sdf begin
            !occursin("-spr", :method)
        end
        m = minimum(sdf1.length)
        beam_best = true
        beam_worst = true
        for r in eachrow(sdf1)
            if r.length == m 
                cons[r.method] += 1
                if !contains(r.method, r"beam") # If this ever occurs, beam was not better than all FastME methods
                    beam_best = false
                else # If this occurs, beam found the best solution
                    beam_worst = false
                end
            end
        end
        if beam_best
            n_beam_was_better += 1
        elseif beam_worst
            n_beam_was_worse += 1
        end
        ## Same with SPR local search
        sdf2 = @rsubset sdf begin
            occursin("-spr", :method)
        end
        m = minimum(sdf2.length)
        beam_best = true
        beam_worst = true
        for r in eachrow(sdf2)
            if r.length == m 
                cons_spr[r.method] += 1
                if !occursin("beam", r.method) # If this ever occurs, beam was not better than all FastME methods
                    beam_best = false
                else # If this occurs, beam found the best solution
                    beam_worst = false
                end
            end
        end
        if beam_best
            n_beam_spr_was_better += 1
        elseif beam_worst
            n_beam_spr_was_worse += 1
        end
    end
    cons_df = DataFrame(method = [], opt = [], opt_perc = [])
    for (m,n) in cons
        push!(cons_df, [m, n, Int(round(n/length(gdf)*100))])
    end
    sort!(cons_df, :opt, rev = true)
    display(cons_df)
    println("Beam was strictly better for ", n_beam_was_better, " datasets (", n_beam_was_better/n_data, ") and strictly worse ", n_beam_was_worse, " times (", n_beam_was_worse/n_data, ")")
    cons_spr_df = DataFrame(method = [], opt = [], opt_perc = [])
    for (m,n) in cons_spr
        push!(cons_spr_df, [m, n, Int(round(n/length(gdf)*100))])
    end
    sort!(cons_spr_df, :opt, rev = true)
    display(cons_spr_df)
    println("Beam+SPR was strictly better for ", n_beam_spr_was_better, " datasets (", n_beam_spr_was_better/n_data, ") and strictly worse ", n_beam_spr_was_worse, " times (", n_beam_spr_was_worse/n_data, ")")
    open("beam_tex_tables.txt", "w") do f
        pretty_table(f, cons_df, backend = Val(:latex), header = ["Method", "Found best", "Percentage"])
        pretty_table(f, cons_spr_df, backend = Val(:latex), header = ["Method", "Found best", "Percentage"])
    end
end