using ProgressMeter
export tree_length, kLPNJ, LP_heuristic, LNS_matheuristic

function tree_length(plm::Matrix, D)
    sum(D[i,j]*2.0^-plm[i,j] for (i,j) in combinations(1:size(D,1),2), init = 0.)*2
end

function nj_select(candidates, τ)
    i, j = first(sort(map(npair -> (NJ_criterion(τ,candidates,npair), npair), combinations(candidates,2)), by = first))[2]
end

function contraction_select(candidates, τ)
    i, j = argmin(tup -> getindex(τ, tup...), combinations(candidates,2))
end

function LP_heuristic(models, _g, D, c, τ; nj_criterion, fastme_postproc = false)
    g = copy(_g)
    sububtnodes = sort(neighbors(_g, c))
    K = length(sububtnodes)
    push!(sububtnodes, c)
    K = degree(g, c)
    while K > 3
        star_tips = Dict(i => n for (i,n) in enumerate(neighbors(g, c)))
        D_check = [D[star_tips[i],star_tips[j]] for i in 1:K, j in 1:K]
        model = models[K]
        set_distances!(model, D_check)
        optimize!(model)
        τ_check = value.(model[:τ]).data
        i_model, j_model = nj_criterion ? nj_select(1:K,τ_check) : contraction_select(1:K,τ_check)
        i, j = star_tips[i_model], star_tips[j_model]
        v = join_neighbors!(g, c, i, j)
        update_distances!(D, g, c, v, i, j, method = :average) # :average is always optimal
        push!(sububtnodes, v)
        K -= 1
    end
    #FastME post-processing
    if fastme_postproc
        g_sub, node_map = induced_subgraph(g, sububtnodes)
        @assert node_map == sububtnodes
        τ_check = path_length_matrix(g_sub)
        D_check = D[node_map[1:K], node_map[1:K]]
        tl_mh = tree_length(τ_check, D_check)
        tmp_D = "tmpLP.txt"
        tmp_tree = "tmpLP.nwk"
        τ_check_FastME = try 
            fastme_local_search(D_check, τ_check, tmp_D, tmp_tree) |> path_length_matrix
        finally
            rm(tmp_D)
            rm(tmp_tree)
            rm(tmp_D*"_fastme_stat.txt")
        end
        tl_fm = tree_length(τ_check_FastME, D_check)
        if tl_fm < tl_mh
            print("(FastME post-processing) ")
        end
        for (i,j) in combinations(1:K,2) # copy to complete τ at correct place
            τ[node_map[i],node_map[j]] = τ[node_map[j],node_map[i]] = τ_check_FastME[i,j]
        end
        while degree(_g,c) > 3
            i,j = contraction_select(neighbors(_g,c), τ)
            v = join_neighbors!(_g, c, i, j)
            update_distances!(τ, _g, c, v, i, j, method = :contraction)
        end
        return _g
    else
        return g
    end
end

function LNS_matheuristic(path, K, tree_path = path*"_tree.nwk"; max_it=Inf, max_time = 60, inittree=false, nj_criterion = false, relax = true, complete = false, fastme_postproc = false)
    time_start = time() 
    time_lim = time_start + max_time
    D = read_distance_matrix(path)
    n = size(D,1)
    τ = Matrix{Float64}(undef, 2n-2, 2n-2) #preallocate once for inplace mutation
    if !inittree
        current_ubt = fastme_local_search(path, tree_path)
        open(tree_path, "w") do f
            write(f, ubt_to_nwk(current_ubt))
        end
    else
        current_ubt = ubtgraph_from_nwk(tree_path)
    end
    PLM = path_length_matrix(current_ubt)
    best_val = tree_length(PLM, D)
    @assert size(D,1) == (nv(current_ubt)+2)÷2
    println("FastME init length: ", best_val)
    it = 0
    last_improve_it = 0
    best_tree_path = tree_path
    models = Dict(k => BMEP_MILP(zeros(k,k); complete) for k in 3:K) # prebuild models of each size, once.
    for model in values(models)
        set_silent(model)
    end
    prog = ProgressUnknown(showspeed = true)
    while it < max_it && time() < time_lim
        it += 1
        next!(prog)
        subtree = sample_neighborhood(current_ubt, K)
        subtree_path_length_matrix!(τ, current_ubt, subtree)
        g_collasped, D_collasped, star_center, node_map = collaspe_to_inner_star(current_ubt, subtree, D)
        new_ubt = LP_heuristic(models, g_collasped, D_collasped, star_center, τ; nj_criterion, fastme_postproc)
        PLM_new = path_length_matrix(new_ubt)
        val = tree_length(PLM_new, D)
        if val < best_val
            best_val = val
            last_improve_it = it
            current_ubt = new_ubt
            time_elapsed = time() - time_start 
            best_tree_path = tree_path*"_$(it)_$time_elapsed"
            open(best_tree_path, "w") do f
                write(f, ubt_to_nwk(current_ubt))
            end
            println(" => $val")
        end
    end
    finish!(prog)
    return best_tree_path
end
