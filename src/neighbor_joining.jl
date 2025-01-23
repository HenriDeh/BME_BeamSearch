using Graphs, Combinatorics, DataStructures
export neighbor_joining, neighbor_joining!, extend_distance, iterated_NJ

function extend_distance(D::Matrix)
    S = similar(D,(size(D).*2 .-2)...)
    S[axes(D)...] .= D
    return S
end

function neighbor_joining(D::AbstractMatrix)
    n = size(D,1)
    g = Graph(n+1)
    c = n+1
    for i in 1:c-1
        add_edge!(g, i, c)
    end
    neighbor_joining!(extend_distance(D),g,c)
end

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

# compute the average distance from v to all neighbors of c, where v is the node that agregates i and j
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
