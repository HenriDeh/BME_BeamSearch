using ProgressMeter
export tree_length, kLPNJ, LP_heuristic, LNS_matheuristic

function tree_length(plm::Matrix, D)
    sum(D[i,j]*2.0^-plm[i,j] for (i,j) in combinations(1:size(D,1),2), init = 0.)*2
end

function nj_select(g, c, τ)
    i, j = first(sort(map(npair -> (NJ_criterion(τ,g,c,npair), npair), combinations(neighbors(g, c),2)), by = first))[2]
end

function contraction_select(g, c, τ)
    i, j = argmin(tup -> getindex(τ, tup...), combinations(neighbors(g,c),2))
end

function LP_heuristic(_g, D, c, τ = similar(D); nj_criterion, iterate = 1, relax = true, complete = false)
    g = copy(_g)
    if complete && !relax
        iterate = nv(g)
    end
    while degree(g, c) > 3
        model = BMEP_MILP(g, D, c, relax = relax, complete = complete)
        set_silent(model)
        optimize!(model) 
        for ((i,j), dist) in value.(model[:τ]).data 
            τ[i,j] = dist
        end
        for _ in 1:iterate
            degree(g, c) > 3 || break
            (i,j) = !nj_criterion ? contraction_select(g,c,τ) : nj_select(g,c,τ)
            v = join_neighbors!(g, c, i, j)
            update_distances!(D, g, c, v, i, j, method = nj_criterion ? :nj : :average)
            update_distances!(τ, g, c, v, i, j, method = nj_criterion ? :nj : :contraction)
        end
    end
    return g
end

function LNS_matheuristic(path, K, tree_path = path*"_tree.nwk"; max_it=Inf, max_time = 60, inittree=false, nj_criterion = false, repair_iterate = 1, relax = true, complete = false)
    time_start = time() 
    time_lim = time_start + max_time
    prog = ProgressUnknown()
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
    while it < max_it && time() < time_lim
        it += 1
        next!(prog)
        subtree = sample_neighborhood(current_ubt, K)
        g_collasped, D_collasped, star_center, node_map = collaspe_to_inner_star(current_ubt, subtree, D)

        new_ubt = LP_heuristic(g_collasped, D_collasped, star_center, τ, nj_criterion=nj_criterion, iterate = repair_iterate, relax = relax, complete = complete)
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
