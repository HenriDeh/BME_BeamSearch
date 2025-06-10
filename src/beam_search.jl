export beam_search, beam_search_lp
using DataStructures

"""
    State

Struct representing a state in the beam search, including the graph, distance matrix, tree length, star tips, and splits.
"""
struct State 
    graph::SimpleGraph
    D::Matrix{Float64}
    l::Float64
    star_tips::Dict{Int,Int}
    splits::Dict{Edge{Int}, Tuple{Set{Int},Set{Int}}}
end

"""
    Cherry

Struct representing a cherry (pair of leaves) in the beam search, with its score and parent state.
"""
struct Cherry
    cherry::Tuple{Int,Int}
    Q::Float64
    parent::State
    splits::Dict{Edge{Int}, Tuple{Set{Int},Set{Int}}}
end

Base.isless(c::Cherry, d::Cherry) = c.Q < d.Q

"""
    reduce_matrix(D, i, j, v, star_tips)

Reduce the distance matrix `D` by merging leaves `i` and `j` into node `v`.
Update the star tips accordingly.
"""
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
    maxmin(x, y)

Return a tuple with the maximum and minimum of `x` and `y`.
"""
maxmin(x,y) = x > y ? (x,y) : (y,x)

"""
    merge_splits(i, j, c, state)

Compute the splits of the candidate tree obtained by merging a cherry `i, j` at state, with `c` being the central node of the star.
"""
function merge_splits(i, j, c, state)
    v = nv(state.graph) + 1
    parent_splits = state.splits
    splitsA1, splitsA2 = parent_splits[Edge(minmax(i, c))]
    splitsB1, splitsB2 = parent_splits[Edge(minmax(j, c))]
    candidate_edges_splits_map = copy(parent_splits)
    new_split = (union(splitsA1, splitsB1), intersect(splitsA2, splitsB2))
    candidate_edges_splits_map[Edge(minmax(c, v))] = new_split
    candidate_edges_splits_map[Edge(maxmin(v, c))] = (new_split[2], new_split[1])
    candidate_edges_splits_map[Edge(minmax(v, i))] = parent_splits[Edge(minmax(i, c))]
    candidate_edges_splits_map[Edge(maxmin(i, v))] = parent_splits[Edge(maxmin(i, c))]
    candidate_edges_splits_map[Edge(minmax(v, j))] = parent_splits[Edge(minmax(j, c))]
    candidate_edges_splits_map[Edge(maxmin(j, v))] = parent_splits[Edge(maxmin(j, c))]
    pop!(candidate_edges_splits_map, Edge(minmax(i, c)))
    pop!(candidate_edges_splits_map, Edge(maxmin(c, i)))
    pop!(candidate_edges_splits_map, Edge(minmax(j, c)))
    pop!(candidate_edges_splits_map, Edge(maxmin(c, j)))
    return candidate_edges_splits_map
end

"""
    beam_search_lp(D_start, B; triangular = true, buneman = false, scale = true, max_coi_cuts = 0, max_length = true, spr = true)

Beam search using linear programming for the Balanced Minimum Evolution problem.
Returns a list of trees.

# Arguments
- `D_start`: Initial distance matrix.
- `B`: Beam width (number of states to keep at each step).
- `triangular`: Use triangle inequalities in the LP.
- `buneman`: Use Buneman constraints in the LP.
- `scale`: Scale objective for numerical stability.
- `max_coi_cuts`: Maximum number of circular ordering inequalities to add.
- `max_length`: Use maximum path length constraint.
- `spr`: Apply FastME local search to each candidate if true.
"""
function beam_search_lp(D_start, B; triangular = true, buneman = false, scale = true, max_coi_cuts = 0, max_length = true, spr = true)
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
    splits_sets = Set{Set{Tuple{Set{Int64}, Set{Int64}}}}()
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
                if length(heap) < B || Q < maximum(heap).Q # If multiple candidates are equivalent but the heap is full, one candidate among the equally worst is arbitrarily discarded to keep the size to B.
                    candidate_edges_splits_map = merge_splits(state.star_tips[i], state.star_tips[j], c, state)
                    candidate_splits = Set(tup for tup in values(candidate_edges_splits_map))
                    if candidate_splits ∉ splits_sets 
                        cherry = Cherry((i,j), Q, state, candidate_edges_splits_map)
                        push!(splits_sets, candidate_splits)
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
            candidate = popmin!(heap)
            # println("The length of the candidate is $(candidate.Q)")
            graph = copy(candidate.parent.graph)
            star_tips = candidate.parent.star_tips
            D = candidate.parent.D
            i_model, j_model = candidate.cherry
            i, j = star_tips[i_model], star_tips[j_model]
            v = join_neighbors!(graph, c, i, j)
            D_check, star_tips = reduce_matrix(copy(D), i_model, j_model, v, star_tips)
            push!(states, State(graph, D_check, candidate.parent.l + D[i_model,j_model]/2, star_tips, candidate.splits))
        end
        K -= 1
    end
    if spr
        return [fastme_local_search(D_start, s.graph) for s in states]
    else
        return [s.graph for s in states]
    end
end

"""
    beam_search(D_start, B; spr = true)

Beam search for the Balanced Minimum Evolution problem using a greedy approach.
Returns a list of trees.

# Arguments
- `D_start`: Initial distance matrix.
- `B`: Beam width (number of states to keep at each step).
- `spr`: Apply FastME local search to each candidate if true.
"""
function beam_search(D_start, B; spr = true)
    _g, c = star_graph(D_start)
    n = only(unique(size(D_start)))
    states = [State(copy(_g), D_start, 0.0, Dict(i => n for (i,n) in enumerate(neighbors(_g, c))), Dict(Edge(i,c) => (Set([i]), setdiff(Set(collect(1:n)), i)) for i in 1:n))]
    for i in 1:n
        states[1].splits[Edge(c,i)] = states[1].splits[Edge(i,c)]
    end
    K = n
    heap = BinaryMinMaxHeap{Cherry}()
    splits_sets = Set{Set{Tuple{Set{Int64}, Set{Int64}}}}()
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
                if length(heap) < B || Q < maximum(heap).Q
                    candidate_edges_splits_map = merge_splits(state.star_tips[i], state.star_tips[j], c, state)
                    candidate_splits = Set(tup for tup in values(candidate_edges_splits_map))
                    if candidate_splits ∉ splits_sets 
                        cherry = Cherry((i,j), Q, state, candidate_edges_splits_map)
                        push!(splits_sets, candidate_splits) 
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
            candidate = popmin!(heap)
            graph = copy(candidate.parent.graph)
            star_tips = candidate.parent.star_tips
            D = candidate.parent.D
            i_model, j_model = candidate.cherry
            i, j = star_tips[i_model], star_tips[j_model]
            v = join_neighbors!(graph, c, i, j)
            D_check, star_tips = reduce_matrix(copy(D), i_model, j_model, v, star_tips)
            push!(states, State(graph, D_check, candidate.parent.l + D[i_model,j_model]/2, star_tips, candidate.splits))
        end
        K -= 1
    end
    if spr
        return [fastme_local_search(D_start, s.graph) for s in states]
    else
        return [s.graph for s in states]
    end
end
