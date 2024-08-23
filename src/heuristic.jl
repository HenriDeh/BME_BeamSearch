using ProgressMeter
export tree_length, kLPNJ, LP_heuristic, LNS_matheuristic

function tree_length(plm::Matrix, D)
    sum(D[i,j]*2.0^-plm[i,j] for (i,j) in combinations(1:size(plm,1),2), init = 0.)*2
end

function NJ_criterion(D, candidates, pair)
    i,j = pair
    (length(candidates)-2)D[i,j] - sum(D[i,k]+D[k,j] for k in candidates if k ∉ (i,j))
end

function nj_select(candidates, τ)
    argmin(pair -> NJ_criterion(τ, candidates, pair), combinations(candidates, 2))
end

function contraction_select(candidates, τ)
    argmin(pair -> getindex(τ, pair...), combinations(candidates,2))
end

# function expNJ_criterion(τ, D, candidates, pair)
#     i,j = pair
#     (length(candidates)-2)D[i,j]*2^-τ[i,j] - sum(D[i,k]*2^-τ[i,k]+D[k,j]*2^-τ[k,j] for k in candidates if k ∉ (i,j))
# end

function expNJ_criterion(τ, D, candidates, pair, taumod = 1)
    i,j = pair
    (length(candidates)-2)D[i,j]*-exp2(-τ[i,j]+taumod) - sum(D[i,k]*-exp2(-τ[i,k]+taumod)+D[k,j]*-exp2(-τ[k,j]+taumod) for k in candidates if k ∉ (i,j))
end

function logNJ_criterion(τ, D, candidates, pair, taumod = 1)
    i,j = pair
    (length(candidates)-2)D[i,j]*log2(τ[i,j]-taumod) - sum(D[i,k]*log2(τ[i,k]-taumod)+D[k,j]*log2(τ[k,j]-taumod) for k in candidates if k ∉ (i,j))
end

function explogNJ_criterion(τ, D, candidates, pair, taumod = 1)
    i,j = pair
    (length(candidates)-2)exp2(D[i,j])*log2(τ[i,j]-taumod) - sum(exp2(D[i,k])*log2(τ[i,k]-taumod)+exp2(D[k,j])*log2(τ[k,j]-taumod) for k in candidates if k ∉ (i,j))
end

function expidNJ_criterion(τ, D, candidates, pair, taumod = 1)
    i,j = pair
    (length(candidates)-2)exp2(D[i,j])*(τ[i,j]-taumod) - sum(exp2(D[i,k])*(τ[i,k]-taumod)+exp2(D[k,j])*(τ[k,j]-taumod) for k in candidates if k ∉ (i,j))
end

function fraclogNJ_criterion(τ, D, candidates, pair, taumod = 1)
    i,j = pair
    (length(candidates)-2)log2(τ[i,j]-taumod) - sum((D[i,k]/D[i,j])*log2(τ[i,k]-taumod)+(D[k,j]/D[i,j])*log2(τ[k,j]-taumod) for k in candidates if k ∉ (i,j))
end

function fracidNJ_criterion(τ, D, candidates, pair, taumod = 1)
    i,j = pair
    (length(candidates)-2)*(τ[i,j]-taumod) - sum((D[i,k]/D[i,j])*(τ[i,k]-taumod)+(D[k,j]/D[i,j])*(τ[k,j]-taumod) for k in candidates if k ∉ (i,j))
end

function expfraclogNJ_criterion(τ, D, candidates, pair, taumod = 1)
    i,j = pair
    (length(candidates)-2)*2*log2(τ[i,j]-taumod) - sum(exp2(D[i,k]/D[i,j])*log2(τ[i,k]-taumod)+exp2(D[k,j]/D[i,j])*log2(τ[k,j]-taumod) for k in candidates if k ∉ (i,j))
end

function expexpNJ_criterion(τ, D, candidates, pair, taumod = 1)
    i,j = pair
    (length(candidates)-2)exp2(D[i,j])*-exp2(-τ[i,j]+taumod) - sum(exp2(D[i,k])*-exp2(-τ[i,k]+taumod)+exp2(D[k,j])*-exp2(-τ[k,j]+taumod) for k in candidates if k ∉ (i,j))
end

function idNJ_criterion(τ, D, candidates, pair, taumod = 1)
    i,j = pair
    (length(candidates)-2)D[i,j]*(τ[i,j]-taumod) - sum(D[i,k]*(τ[i,k]-taumod)+D[k,j]*(τ[k,j]-taumod) for k in candidates if k ∉ (i,j))
end

function tau_select(candidates, τ, D)
    argmin(pair -> fracidNJ_criterion(τ, D, candidates, pair, 1.0), combinations(candidates, 2))
end

function tau_select(candidates, τ)
    n = length(candidates)
    tau_select(candidates, τ, ones(n,n))
end

function LP_heuristic(_g, D, c; nj_criterion, fastme_postproc = false, complete = true, scale = false) # Stand alone Method
    n = length(neighbors(_g, c))
    models = Dict(k => BMEP_MILP(zeros(k,k); complete, scale) for k in 3:n)
    for model in models
        set_silent(model.second)
    end
    LP_heuristic(models, _g, D, c, similar(D); nj_criterion, fastme_postproc, scale)
end

function LP_heuristic(models, _g, D, c, τ; nj_criterion, fastme_postproc = false, scale = false, max_length = true)
    g = copy(_g)
    sububtnodes = sort(neighbors(_g, c))
    push!(sububtnodes, c)
    K = degree(g, c)
    while K > 3
        star_tips = Dict(i => n for (i,n) in enumerate(neighbors(g, c)))
        D_check = [D[star_tips[i],star_tips[j]] for i in 1:K, j in 1:K]
        model = models[K]
        set_distances!(model, D_check; scale)
        optimize!(model)
        candidates = 1:K
        tau = [i == j ? Inf : value.(model[:τ][i,j]) for i in candidates, j in candidates]
        # x = value.(model[:x])
        # tau = [i == j ? 0. : argmax(l -> x[i, j, l], 2:(max_length ? min(K-1, ceil(log2(K-1)^2)) : K-1)) for i in 1:K, j in 1:K]
        # if nj_criterion
        #     tauheap = sort([pair => logNJ_criterion(tau, D_check,  candidates, pair, .0) for pair in combinations(candidates, 2)], by = p -> p.second)
        #     i_model, j_model = first(sort([pair => NJ_criterion(D, candidates, pair) for (pair, d) in tauheap[1:3]], by= p->p.second)).first
        # else
        #     tauheap =sort([pair => tau[pair...] for pair in combinations(candidates, 2)], by = p -> p.second)
        #     i_model, j_model = first(tauheap).first
        # end
        tauheap = sort([pair => tau[pair...] for pair in combinations(candidates, 2)], by = p -> p.second)
        k = max(6, findfirst(p -> p.second > 2, tauheap) - 1)
        i_model, j_model = first(sort([pair => logNJ_criterion(tau, D_check,  candidates, pair, 1.0) for (pair, t) in tauheap[1:k]], by = p -> p.second)).first
         #first(sort([pair => NJ_criterion(D, candidates, pair) for (pair, d) in taunj[1:3]], by= p->p.second)).first

        # cbests = Set{Vector{Int}}()
        # njbests = Set{Vector{Int}}()
        
        # i_model, j_model = nj_criterion ? tau_select(1:K, tau, D_check) : contraction_select(1:K,tau)
        i, j = star_tips[i_model], star_tips[j_model]
        v = join_neighbors!(g, c, i, j)
        update_distances!(D, g, c, v, i, j, method = :average) # :average is always optimal
        push!(sububtnodes, v)
        K -= 1
    end
    #FastME post-processing
    if fastme_postproc
        n = length(leaves(g))
        # K = degree(g, c)
        # g_sub, node_map = induced_subgraph(g, sububtnodes)
        # @assert node_map == sububtnodes
        τ_check = path_length_matrix(g)
        # D_check = D[node_map[1:n], node_map[1:n]]
        D_check = D[1:n,1:n]
        # tl_mh = tree_length(τ_check, D_check)
        tmp_D = "tmpLP.txt"
        tmp_tree = "tmpLP.nwk"
        # τ_check_FastME = try 
        #     fastme_local_search(D_check, τ_check, tmp_D, tmp_tree) |> path_length_matrix
        # finally
        #     rm(tmp_D)
        #     rm(tmp_tree)
        #     rm(tmp_D*"_fastme_stat.txt")
        # end
        # tl_fm = tree_length(τ_check_FastME, D_check)
        # if tl_fm < tl_mh
        #     print(" (FME post-proc: $tl_mh => $tl_fm) ")
        # end
        # for (i,j) in combinations(1:K,2) # copy to complete τ at correct place
        #     τ[node_map[i],node_map[j]] = τ[node_map[j],node_map[i]] = τ_check_FastME[i,j]
        # end
        # while degree(_g,c) > 3
        #     i,j = contraction_select(neighbors(_g,c), τ)
        #     v = join_neighbors!(_g, c, i, j)
        #     update_distances!(τ, _g, c, v, i, j, method = :contraction)
        # end
        g = try
            fastme_local_search(D_check, τ_check, tmp_D, tmp_tree)
        finally
            rm(tmp_D)
            rm(tmp_tree)
            rm(tmp_D*"_fastme_stat.txt")
        end
    end
    return g
end

function LNS_matheuristic(path, K, tree_path = path*"_tree.nwk"; max_it=Inf, max_time = 60, inittree=false, nj_criterion = false, relax = true, complete = false, fastme_postproc = false, matheuristic = true)
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
    print("Building models...")
    models = Dict(k => BMEP_MILP(zeros(k,k); complete, relax) for k in 3:(matheuristic ? K : 0)) # prebuild models of each size, once.
    for model in values(models)
        set_silent(model)
    end
    println("\b\b\b done.")
    prog = ProgressUnknown(showspeed = true)
    while it < max_it && time() < time_lim
        it += 1
        next!(prog)
        subtree = sample_neighborhood(current_ubt, K)
        g_collasped, D_collasped, star_center, node_map = collaspe_to_inner_star(current_ubt, subtree, D)
        if matheuristic
            new_ubt = LP_heuristic(models, g_collasped, D_collasped, star_center, τ; nj_criterion, fastme_postproc)
        else
            new_ubt = FastME_heuristic(path, g_collasped, D_collasped, star_center, τ)
        end
        val = tree_length(path_length_matrix(new_ubt), D)
        if val < best_val
            best_val = val
            last_improve_it = it
            current_ubt = new_ubt
            time_elapsed = round(time() - time_start)
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

function FastME_heuristic(path, _g, D, c, τ)
    g_best = copy(_g)
    star_tips = Dict(i => n for (i,n) in enumerate(neighbors(g_best, c)))
    K = degree(g_best, c)
    D_check = [i == j ? 0. : D[minmax(star_tips[i], star_tips[j])...] for i in 1:K, j in 1:K]
    tmp_D = "tmpLP.txt"
    tmp_tree = "tmpLP.nwk"
    best_val = Inf
    try 
        for method in ["b","i","o","n","u"]
            g = copy(_g)
            τ_check_FastME = try
                fastme_local_search(D_check, tmp_D, tmp_tree, methods = [method], verbose = false, spr = false) |> path_length_matrix
            catch e
                continue
            end
            for (i,j) in combinations(1:K,2) # copy to complete τ at correct place
                τ[star_tips[i],star_tips[j]] = τ[star_tips[j],star_tips[i]] = τ_check_FastME[i,j]
            end
            while degree(g,c) > 3
                i,j = contraction_select(neighbors(g,c), τ)
                v = join_neighbors!(g, c, i, j)
                update_distances!(τ, g, c, v, i, j, method = :contraction)
            end
            g = try 
                open("tmpLP.nwk", "w") do f
                    write(f, ubt_to_nwk(g))
                end
                fastme_local_search(path, "tmpLP.nwk", inittree = true)
            finally
                rm("tmpLP.nwk")
            end
            val = tree_length(path_length_matrix(g), D)
            if val < best_val
                g_best = g
                best_val = val
            end
        end    
    finally
        rm(tmp_D)
        rm(tmp_tree, force = true)
        rm(tmp_D*"_fastme_stat.txt")
    end
    return g_best
end

