export beam_search, beam_search_lp
using DataStructures

struct State 
    graph::SimpleGraph
    D::Matrix{Float64}
    l::Float64
    star_tips::Dict{Int,Int}
    splits::Dict{Edge{Int}, Tuple{Set{Int},Set{Int}}}
end

struct Cherry
    cherry::Tuple{Int,Int}
    Q::Float64
    parent::State
    splits::Dict{Edge{Int}, Tuple{Set{Int},Set{Int}}}
end

Base.isless(c::Cherry, d::Cherry) = c.Q < d.Q

function reduce_matrix(D, i, j, v, star_tips)
    n = only(unique(size(D)))
    i, j = minmax(i,j)
    D[i,:] .= D[:,i] .= (D[i,:] .+ D[j,:]) ./2
    D_check = reshape([a == b ? 0.0 : D[a,b] for a in 1:n, b in 1:n if a != j && b != j], n-1,n-1)
    node_map = [k < j ? k : k + 1 for k in 1:n-1] # new star_tips
    star_tips = Dict(k => k == i ? v : star_tips[node_map[k]] for k in 1:n-1)
    return D_check, star_tips
end

"""
    cherry_splits(cherry, state, c)

    Computes the splits of the child tree obtained by merging a cherry.
"""
function cherry_splits(i, j, c, state)
    v = nv(state.graph) + 1
    parent_splits = state.splits
    splitsA1, splitsA2 = parent_splits[Edge(minmax(i, c))]
    splitsB1, splitsB2 = parent_splits[Edge(minmax(j, c))]
    child_splits = copy(parent_splits)
    child_splits[Edge(minmax(c, v))] = (union(splitsA1, splitsB1), intersect(splitsA2, splitsB2))
    child_splits[Edge(minmax(v, i))] = parent_splits[Edge(minmax(i, c))]
    child_splits[Edge(minmax(v, j))] = parent_splits[Edge(minmax(j, c))]
    pop!(child_splits, Edge(minmax(i, c)))
    pop!(child_splits, Edge(minmax(j, c)))
    return child_splits
end


"""
    push_to_heap!(heap, B, cherry, state, c, splits_sets)

    Pushes a cherry to the heap if it is better than the worst cherry in the heap. 
    If the heap is not full, the cherry is pushed to the heap. 
    Otherwise, the cherry is pushed to the heap if it is better than the worst cherry in the heap.
    Pushing only occurs if the splits of the child tree are not already in the heap.

    A merge consist in adding one split to the child tree. (i, c) and (j, c) and (v, i) and (v, j).
    The new split is (v, c) and is computed by merging the i and j sides of the splits of the parent tree.
"""
function push_to_heap!(heap, B, cherry, splits_sets)
   
end

function beam_search_lp(D_start, B; triangular = true, buneman = false, scale = true, max_coi_cuts = 0, max_length = true, spr = true, alltwos = false)
    _g, c = star_graph(D_start)
    n = only(unique(size(D_start)))
    models = Dict(k => BMEP_MILP(zeros(k,k); triangular, buneman, scale, max_length) for k in 3:n)
    for model in models
        set_silent(model.second)
    end
    states = [State(copy(_g), D_start, 0.0, Dict(i => n for (i,n) in enumerate(neighbors(_g, c))), Dict(Edge(i,c) => (Set([i]), setdiff(Set(collect(1:n)), i)) for i in 1:n))]
    for i in 1:n
        states[1].splits[Edge(c,i)] = states[1].splits[Edge(i,c)]
    end
    K = n
    heap = BinaryMinMaxHeap{Cherry}()
    splits_sets = Set([values(only(states).splits)])
    while K > 3
        empty!(heap)
        empty!(splits_sets)
        for state in states
            D = state.D
            model = models[K]
            set_distances!(model, D; scale)
            solve_BME_model!(model; max_coi_cuts)
            LB = value(model[:tree_length])
            # println("The length of the lower bound is $(LB + state.l)")
            τ = value.(model[:τ])
            for (i,j) in combinations(1:K,2)
                lb_delta = D[i,j]/2 - D[i,j]/exp2(τ[i,j]-1)
                Q = LB + state.l + lb_delta
                if isempty(heap) || (Q < maximum(heap).Q || length(heap) < B) || (alltwos && Q == maximum(heap).Q) 
                    child_splits = cherry_splits(state.star_tips[i], state.star_tips[j], c, state) # If multiple childs are equivalent but the heap is full, childs are arbitrarily discarded to keep the size to B.
                    if child_splits ∉ splits_sets 
                        cherry = Cherry((i,j), Q, state, child_splits)
                        push!(splits_sets, values(cherry.splits))
                        push!(heap, cherry)
                        if length(heap) > B 
                            popmax!(heap)
                        end
                    end 
                end 
            end
        end
        empty!(states)
        while length(states) < B && !isempty(heap)
            child = popmin!(heap)
            # println("The length of the child is $(child.Q)")
            gchild = copy(child.parent.graph)
            star_tips = child.parent.star_tips
            D = child.parent.D
            i_model, j_model = child.cherry
            i, j = star_tips[i_model], star_tips[j_model]
            v = join_neighbors!(gchild, c, i, j)
            D_check, star_tips = reduce_matrix(copy(D), i_model, j_model, v, star_tips)
            push!(states, State(gchild, D_check, child.parent.l + D[i_model,j_model]/2, star_tips, child.splits))
        end
        K -= 1
    end
    # println("The length of the solution is $(states[1].l + sum(states[1].D)/4)")
    if spr
        return [fastme_local_search(D_start, s.graph) for s in states]
    else
        return [s.graph for s in states]
    end
end

function beam_search(D_start, B; spr = true)
    _g, c = star_graph(D_start)
    n = only(unique(size(D_start)))
    states = [State(copy(_g), D_start, 0.0, Dict(i => n for (i,n) in enumerate(neighbors(_g, c))), Dict(Edge(i,c) => (Set([i]), setdiff(Set(collect(1:n)), i)) for i in 1:n))]
    for i in 1:n
        states[1].splits[Edge(c,i)] = states[1].splits[Edge(i,c)]
    end
    K = n
    heap = BinaryMinMaxHeap{Cherry}()
    splits_sets = Set([values(only(states).splits)])
    while K > 3
        empty!(heap)
        empty!(splits_sets)
        for state in states
            D = state.D
            col_sums = sum(D, dims = 1)
            sum_D = sum(col_sums)
            for (i,j) in combinations(1:K, 2)
                star_tree_length = (sum_D - col_sums[i] - col_sums[j])/(2(K-2))
                Q = state.l + D[i,j]/2 + star_tree_length
                if isempty(heap) || (Q < maximum(heap).Q || length(heap) < B) 
                    child_splits = cherry_splits(state.star_tips[i], state.star_tips[j], c, state) # If multiple childs are equivalent but the heap is full, childs are arbitrarily discarded to keep the size to B.
                    if child_splits ∉ splits_sets 
                        cherry = Cherry((i,j), Q, state, child_splits)
                        push!(splits_sets, values(cherry.splits))
                        push!(heap, cherry)
                        if length(heap) > B 
                            popmax!(heap)
                        end
                    end
                end 
            end
        end
        empty!(states)
        while length(states) < B && !isempty(heap)
            child = popmin!(heap)
            gchild = copy(child.parent.graph)
            star_tips = child.parent.star_tips
            D = child.parent.D
            i_model, j_model = child.cherry
            i, j = star_tips[i_model], star_tips[j_model]
            v = join_neighbors!(gchild, c, i, j)
            D_check, star_tips = reduce_matrix(copy(D), i_model, j_model, v, star_tips)
            push!(states, State(gchild, D_check, child.parent.l + D[i_model,j_model]/2, star_tips, child.splits))
        end
        K -= 1
    end
    if spr
        return [fastme_local_search(D_start, s.graph) for s in states]
    else
        return [s.graph for s in states]
    end
end
