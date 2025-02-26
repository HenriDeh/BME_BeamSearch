export beam_search, beam_search_lp

function reduce_matrix(D, i, j, v, star_tips)
    n = only(unique(size(D)))
    i, j = minmax(i,j)
    D[i,:] .= D[:,i] .= (D[i,:] .+ D[j,:]) ./2
    D_check = reshape([a == b ? 0.0 : D[a,b] for a in 1:n, b in 1:n if a != j && b != j], n-1,n-1)
    node_map = [k < j ? k : k + 1 for k in 1:n-1] # new star_tips
    star_tips = Dict(k => k == i ? v : star_tips[node_map[k]] for k in 1:n-1)
    return D_check, star_tips
end

function beam_search_lp(D_start, B; triangular = true, buneman = false, scale = true, max_coi_cuts = 0, max_length = true, spr = true)
    _g, c = star_graph(D_start)
    n = only(unique(size(D_start)))
    models = Dict(k => BMEP_MILP(zeros(k,k); triangular, buneman, scale, max_length) for k in 3:n)
    for model in models
        set_silent(model.second)
    end
    states = [(graph = copy(_g), D = D_start, l = 0.0, star_tips = Dict(i => n for (i,n) in enumerate(neighbors(_g, c))))]
    K = n
    while K > 3
        heap = @NamedTuple{cherry::Tuple{Int,Int}, Q::Float64, parent::eltype(states)}[] # the heap could be optimized with an actual heap DS.
        for state in states
            D = state.D
            model = models[K]
            set_distances!(model, D; scale)
            solve_BME_model!(model; max_coi_cuts)
            candidates = 1:K
            LB = value(model[:tree_length])
            τ = value.(model[:τ])
            for (i,j) in combinations(candidates, 2)
                lb_delta = D[i,j]/2 - D[i,j]/exp2(τ[i,j]-1)
                push!(heap, (cherry = (i,j), Q = LB + state.l + lb_delta, parent = state))
            end
        end
        sort!(heap, rev = false, by = p -> p.Q)
        empty!(states)
        split_sets = Set{Set{Tuple{Set{Int}, Set{Int}}}}()
        idx = 0
        while length(states) < B && idx <= length(heap)
            idx += 1
            child = heap[idx]
            gchild = copy(child.parent.graph)
            star_tips = child.parent.star_tips
            D = child.parent.D
            i_model, j_model = child.cherry
            i, j = star_tips[i_model], star_tips[j_model]
            v = join_neighbors!(gchild, c, i, j)
            splits = compute_splits(gchild)
            if splits ∉ split_sets
                push!(split_sets, splits)
                D_check, star_tips = reduce_matrix(copy(D), i_model, j_model, v, star_tips) # averages
                push!(states, (graph = gchild, D = D_check, l = child.parent.l + D[i_model,j_model]/2, star_tips = star_tips))
            end
        end
        K -= 1
    end
    if spr
        return [fastme_local_search(D_start, s.graph) for s in states]
    else
        return [s.graph for s in states]
    end
end

function beam_search(D_start, B; spr = true)
    _g, c = star_graph(D_start)
    n = only(unique(size(D_start)))
    states = [(graph = copy(_g), D = D_start, l = 0.0, star_tips = Dict(i => n for (i,n) in enumerate(neighbors(_g, c))))]
    K = n
    while K > 3
        heap = @NamedTuple{cherry::Tuple{Int,Int}, Q::Float64, parent::eltype(states)}[] # the heap could be optimized with an actual heap DS.
        for state in states
            D = state.D
            candidates = 1:K
            col_sums = sum(D, dims = 1)
            sum_D = sum(col_sums)
            for (i,j) in combinations(candidates, 2)
                star_tree_length = (sum_D - col_sums[i] - col_sums[j])/(2(K-2))
                push!(heap, (cherry = (i,j), Q = state.l + D[i,j]/2 + star_tree_length, parent = state))
            end
        end
        sort!(heap, rev = false, by = p -> p.Q)
        empty!(states)
        split_sets = Set{Set{Tuple{Set{Int}, Set{Int}}}}()
        idx = 0
        while length(states) < B && idx <= length(heap)
            idx += 1
            child = heap[idx]
            gchild = copy(child.parent.graph)
            star_tips = child.parent.star_tips
            D = child.parent.D
            i_model, j_model = child.cherry
            i, j = star_tips[i_model], star_tips[j_model]
            v = join_neighbors!(gchild, c, i, j)
            splits = compute_splits(gchild)
            if splits ∉ split_sets
                push!(split_sets, splits)
                D_check, star_tips = reduce_matrix(copy(D), i_model, j_model, v, star_tips) # averages
                push!(states, (graph = gchild, D = D_check, l = child.parent.l + D[i_model,j_model]/2, star_tips = star_tips))
            end
        end
        K -= 1
    end
    if spr
        return [fastme_local_search(D_start, s.graph) for s in states]
    else
        return [s.graph for s in states]
    end
end
