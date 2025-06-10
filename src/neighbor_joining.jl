using Graphs, Combinatorics, DataStructures
export neighbor_joining, neighbor_joining!, extend_distance, iterated_NJ

"""
    extend_distance(D::Matrix)

Extend the distance matrix `D` to the appropriate size for tree construction.
"""
function extend_distance(D::Matrix)
    S = similar(D,(size(D).*2 .-2)...)
    S[axes(D)...] .= D
    return S
end

"""
    neighbor_joining(D::AbstractMatrix)

Construct a tree using the neighbor joining algorithm from distance matrix `D`.
"""
function neighbor_joining(D::AbstractMatrix)
    n = size(D,1)
    g = Graph(n+1)
    c = n+1
    for i in 1:c-1
        add_edge!(g, i, c)
    end
    neighbor_joining!(extend_distance(D),g,c)
end

"""
    neighbor_joining!(D, g, c)

In-place neighbor joining on graph `g` with distance matrix `D` and center node `c`.

# Arguments
- `D`: Distance matrix (should be square and large enough for all nodes in `g`).
- `g`: Graph object, initially a star with leaves and center `c`.
- `c`: The center node of the star (integer node index).
"""
function neighbor_joining!(D,g,c) # D is a distance matrix (2n-2 x 2n-2), g is a graph containing a star with n branches and internal node c = n+1
    nj_crit = BinaryHeap(Base.By(last), [(maximum(npair),minimum(npair)) => NJ_criterion(D,neighbors(g, c),npair) for npair in combinations(neighbors(g, c),2)])
    joined = Set{Int}()
    while degree(g, c) > 3
        #println(degree(g,c))
        (i,j), dist = pop!(nj_crit)
        if i ∉ joined && j ∉ joined
            v = join_neighbors!(g, c, i, j)
            push!(joined, i)
            push!(joined, j)
            update_distances!(D, g, c, v, i, j, method = :average)
            for k in neighbors(g,c) 
                if k ≠ v 
                    push!(nj_crit, (v,k) => NJ_criterion(D,neighbors(g,c),(v,k)))
                end
            end
        end
    end
    return g
end

"""
    join_neighbors!(g, c, i, j)

Join neighbors `i` and `j` at center `c` in graph `g`, returning the new node.
"""
function join_neighbors!(g, c, i, j)
    rem_edge!(g, c, i)
    rem_edge!(g, c, j)
    add_vertex!(g)
    v = nv(g)
    add_edge!(g, v, i)
    add_edge!(g, v, j)
    add_edge!(g, v, c)
    return v
end

"""
    update_distances!(D, g, c, v, i, j; method = :average)

Update the distance matrix `D` after joining nodes `i` and `j` into `v` at center `c`.

# Arguments
- `D`: Distance matrix to update (modified in-place).
- `g`: Graph object.
- `c`: Center node of the star.
- `v`: New node created by joining `i` and `j`.
- `i`, `j`: Nodes being joined.
- `method`: Method for updating distances (`:average`, `:dij`, or `:contraction`).
"""
function update_distances!(D, g, c, v, i, j; method = :average)
    if method == :dij
        reduction = D[i, j]
    elseif method == :average
        reduction = 0.
    elseif method == :contraction
        reduction = 2.
    else
        throw(ArgumentError("Unknown reduction method $(string(method))"))
    end    
    for k in neighbors(g, c) 
        if k ≠ v
            D[k,v] = D[v,k] = (D[i,k] + D[k,j] - reduction)/2
        end
    end
    D[v,v] = 0.
end
