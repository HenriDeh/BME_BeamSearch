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

function LP_heuristic(_g, D, c, τ = similar(D); nj_criterion, iterate = 1, relax = true, exact = false)
    g = copy(_g)
    if exact && !relax
        iterate = nv(g)
    end
    while degree(g, c) > 3
        if exact
            model = MIP_complete(g, D, c, relax = relax)
        else
            model = MIP_reduced(g, D, c, relax = relax)
        end
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

function LNS_matheuristic(path, K, tree_path = path*"_tree.nwk"; max_it=Inf, max_time = 60, fastme_it = Inf, inittree=false, nj_criterion = false, repair_iterate = 1, relax = true, exact = false, nj_repair = false, check_radius = false)
    time_lim = time() + max_time
    prog = ProgressUnknown()
    D = read_distance_matrix(path)
    n = size(D,1)
    τ = Matrix{Float64}(undef, 2n-2, 2n-2) #preallocate once for inplace mutation
    if !inittree
        current_ubt = fastme_local_search(path, tree_path)
        #current_ubt = neighbor_joining(D)
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
    # Lower bound
    if check_radius
        print("Computing lowerbound: ")
        gmh = Graph(n+1)
        c = n+1
        for i in 1:c-1
            add_edge!(gmh, i, c)
        end
        D_ = extend_distance(D);
        τmh = similar(D_);
        tau_tilde = similar(D)
        while degree(gmh, c) > 3
            #model = MIP_complete(gmh, D_, c, relax = true)
            model = MIP_reduced(gmh, D_, c, relax = true)
            set_silent(model)
            optimize!(model) 
            for ((i,j), dist) in value.(model[:τ]).data 
                τmh[i,j] = dist
                τmh[i,i] = 0.
                τmh[j,j] = 0.
            end
            if degree(gmh, c) == n
                tau_tilde .= τmh[1:n,1:n]
            end
            (i,j) = !nj_criterion ? contraction_select(gmh,c,τmh) : nj_select(gmh,c,τmh)
            v = join_neighbors!(gmh, c, i, j)
            update_distances!(D_, gmh, c, v, i, j, method = nj_criterion ? :nj : :average)
            update_distances!(τmh, gmh, c, v, i, j, method = nj_criterion ? :nj : :contraction)
        end
        lb = tree_length(tau_tilde, D)
        println(lb," ", extrema(tau_tilde))
        PLM_MH = path_length_matrix(gmh)
        tl_mh = tree_length(PLM_MH, D)
        print("Matheuristic found initial tree with length $tl_mh, ")
        if tl_mh < best_val
            current_ubt = gmh
            best_val = tl_mh
            PLM = PLM_MH
            open(tree_path, "w") do f
                write(f, ubt_to_nwk(current_ubt))
            end
            println("using it as initial tree.")
        else
            println("keeping FastME tree.")
        end
        radius = maximum(abs.(PLM .- tau_tilde[1:n,1:n]))
        println("Radius = $radius")
        if radius <= 0.5
            println("Solution is within safety radius, thus optimal, skipping Large Neighborhood Search")
            return tree_path
        end
        println("LP Gap = ", 1-lb/best_val)
    end
    it = 0
    its_since_improve = 0
    fastme_can_improve = false
    neighborhood = Rand_Neighborhood(current_ubt) 
    last_improve_it = 0
    while it < max_it && time() < time_lim
        it += 1
        next!(prog)
        g_collasped, D_collasped, star_center, node_map = collaspe_to_inner_star(current_ubt, K, D, neighborhood = neighborhood)
        new_ubt = LP_heuristic(g_collasped, D_collasped, star_center, τ, nj_criterion=nj_criterion, iterate = repair_iterate, relax = relax, exact = exact)
        PLM_new = path_length_matrix(new_ubt)
        val = tree_length(PLM_new, D)
        if val < best_val
            best_val = val
            last_improve_it = it
            current_ubt = new_ubt
            open(tree_path*"_$it", "w") do f
                write(f, ubt_to_nwk(current_ubt))
            end
            its_since_improve = 0
            fastme_can_improve = true 
            if check_radius
                radius = maximum(abs.(PLM .- tau_tilde[1:n,1:n]))
                println(" => $val LP Gap = ", lb/best_val-1, " Radius = $radius")
                if radius <= 0.5
                    println("Solution is within safety radius, thus optimal, stopping early")
                    return tree_path*"_$it"
                end
            else
                println(" => $val")
            end
        else
            its_since_improve +=1
            if its_since_improve >= fastme_it && fastme_can_improve
                print(" Performing SPR ")
                new_ubt = fastme_local_search(path, tree_path, inittree = true)
                val = tree_length(path_length_matrix(new_ubt), D)
                if val < best_val
                    println("=> ", val)
                    best_val = val
                    current_ubt = new_ubt  
                else
                    println(" (no improvement found)")
                    open(tree_path, "w") do f
                        write(f, ubt_to_nwk(current_ubt))
                    end
                end
                fastme_can_improve = false # don't retry fastme until the tree improved
            end
        end
    end
    finish!(prog)
    return tree_path * (last_improve_it == 0 ? "" : "_$last_improve_it")
end