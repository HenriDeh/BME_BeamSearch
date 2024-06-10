export Rand_Neighborhood, Weighted_Neighborhood, BFS_Neighborhood, DFS_Neighborhood

struct Rand_Neighborhood
end
Rand_Neighborhood(::Graph) = Rand_Neighborhood()

function (r::Rand_Neighborhood)(g, K)
    start = rand(inner_nodes(g))
    selected = Set([start])
    leaves = Set(filter(n->degree(g,n) == 1, neighbors(g,start)))
    candidates = Set{Int}(filter(n->degree(g,n) == 3, neighbors(g,start)))
    while length(selected) < K-2
        new_node = rand(candidates)
        push!(selected, new_node)
        pop!(candidates, new_node)
        union!(candidates, setdiff(filter(n->degree(g,n) == 3, neighbors(g,new_node)), selected)) #add its neighbors in candidates if they have degree 3
        union!(leaves, filter(n->degree(g,n) == 1, neighbors(g,new_node)))  #neighbors with degree 1 are leaves of g
    end
    selected, leaves, candidates
end

mutable struct Weighted_Neighborhood
    weights::OffsetVector{Float64, Vector{Float64}}
end
function Weighted_Neighborhood(g::Graph) 
    nplus1 = first(inner_nodes(g))
    rinner = nplus1:nv(g)
    Weighted_Neighborhood(OffsetVector(ones(length(rinner)), rinner))
end


function (r::Weighted_Neighborhood)(g, K)
    start = sample(inner_nodes(g), FrequencyWeights(r.weights.parent))
    selected = Set([start])
    leaves = Set(filter(n->degree(g,n) == 1, neighbors(g,start)))
    candidates = Set{Int}(filter(n->degree(g,n) == 3, neighbors(g,start)))
    while length(selected) < K-2
        cs = collect(candidates)
        ws = [r.weights[c] for c in cs]
        new_node = sample(cs, FrequencyWeights(ws))
        push!(selected, new_node)
        pop!(candidates, new_node)
        union!(candidates, setdiff(filter(n->degree(g,n) == 3, neighbors(g,new_node)), selected)) #add its neighbors in candidates if they have degree 3
        union!(leaves, filter(n->degree(g,n) == 1, neighbors(g,new_node)))  #neighbors with degree 1 are leaves of g
    end
    selected, leaves, candidates
end

update_weights!(::Any, node_map) = nothing
function update_weights!(r::Weighted_Neighborhood, node_map, selected)
    new_weights = similar(r.weights)
    new_weights .= 1.
    inner_nodes = eachindex(new_weights)
    inner_node_map = node_map[first(inner_nodes):end]
    offset = new_weights.offset
    for (nm_idx, old_idx) in enumerate(inner_node_map)
        new_weights[nm_idx + offset] = r.weights[old_idx] *10
    end
    r.weights = new_weights
    return nothing
end

mutable struct BFS_Neighborhood
    root::Int
end
BFS_Neighborhood(g::Graph) = BFS_Neighborhood(first(inner_nodes(g)))

function (bfsn::BFS_Neighborhood)(g, K)
    start = bfsn.root
    bfsn.root += 1
    if bfsn.root > nv(g)
        bfsn.root = first(inner_nodes(g))
    end
    selected = Set([start])
    leaves = Set(filter(n->degree(g,n) == 1, neighbors(g,start)))
    candidates = sort(filter(n->degree(g,n) == 3, neighbors(g,start)))
    while length(selected) < K-2
        new_node = popfirst!(candidates)
        push!(selected, new_node)
        for nei = neighbors(g,new_node)
            if degree(g,nei) == 1
                push!(leaves, nei)
            elseif nei ∉ selected
                push!(candidates, nei)
            end
        end
    end
    return selected, leaves, candidates
end

mutable struct DFS_Neighborhood
    root::Int
end
DFS_Neighborhood(g::Graph) = DFS_Neighborhood(first(inner_nodes(g)))

function (dfsn::DFS_Neighborhood)(g, K)
    start = dfsn.root
    dfsn.root += 1
    if dfsn.root > nv(g)
        dfsn.root = first(inner_nodes(g))
    end
    selected = Set([start])
    leaves = Set(filter(n->degree(g,n) == 1, neighbors(g,start)))
    candidates = Stack{Int}()
    for c in sort(filter(n->degree(g,n) == 3, neighbors(g,start)))
        push!(candidates,c)
    end

    while length(selected) < K-2
        new_node = pop!(candidates)
        push!(selected, new_node)
        for nei = neighbors(g,new_node)
            if degree(g,nei) == 1
                push!(leaves, nei)
            elseif nei ∉ selected
                push!(candidates, nei)
            end
        end
    end
    return selected, leaves, candidates
end

mutable struct Distance_Neighborhood
    ordering::Vector{Vector{Int}}
    τ::Matrix{Float64}
end

function Distance_Neighborhood(g::Graph, D::Matrix, τ)
    ordering = sort(collect(combinations(leaves(g), 2)), by = v -> D[v...])
    Distance_Neighborhood(ordering, τ)
end

function (dbn::Distance_Neighborhood)(g, K)
    τmin = nv(g)
    rinner = inner_nodes(g)
    weights = OffsetVector(ones(length(rinner)), rinner)
    fw = floyd_warshall_shortest_paths(g)
    τ = dbn.τ
    for (i,j) in Iterators.reverse(dbn.ordering)
        if τ[i,j] > τmin
            parent = fw.parents[i,j]
            while parent != i
                weights[parent] += 1
                parent = fw.parents[i, parent]
            end
        else
            τmin = τ[i,j]
        end
    end

    start = sample(inner_nodes(g), FrequencyWeights(weights.parent))
    selected = Set([start])
    leaves = Set(filter(n->degree(g,n) == 1, neighbors(g,start)))
    candidates = Set{Int}(filter(n->degree(g,n) == 3, neighbors(g,start)))
    while length(selected) < K-2
        cs = collect(candidates)
        ws = [weights[c] for c in cs]
        new_node = sample(cs, FrequencyWeights(ws))
        push!(selected, new_node)
        pop!(candidates, new_node)
        union!(candidates, setdiff(filter(n->degree(g,n) == 3, neighbors(g,new_node)), selected)) #add its neighbors in candidates if they have degree 3
        union!(leaves, filter(n->degree(g,n) == 1, neighbors(g,new_node)))  #neighbors with degree 1 are leaves of g
    end
    selected, leaves, candidates
end