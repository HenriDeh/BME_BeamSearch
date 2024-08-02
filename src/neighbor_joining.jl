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
    nj_crit = BinaryHeap(Base.By(last), [(maximum(npair),minimum(npair)) => NJ_criterion(D,g,c,npair) for npair in combinations(neighbors(g, c),2)])
    joined = Set{Int}()
    while degree(g, c) > 3
        #println(degree(g,c))
        (i,j), dist = pop!(nj_crit)
        if i ∉ joined && j ∉ joined
            v = join_neighbors!(g, c, i, j)
            push!(joined, i)
            push!(joined, j)
            update_distances!(D, g, c, v, i, j)
            for k in neighbors(g,c) 
                if k ≠ v 
                    push!(nj_crit, (v,k) => NJ_criterion(D,g,c,(v,k)))
                end
            end
        end
    end
    return g
end

function NJ_criterion(D, candidates, pair) #Q criterion p. 16 Catanzaro et al. 2020
    i,j = pair
    (length(candidates))*distance(D,i,j) - sum(distance(D,i,k)+distance(D,k,j) for k in candidates if k ∉ (i,j))
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
function update_distances!(D, g, c, v, i, j; method = :nj)
    reduction = (method == :nj ? distance(D,i,j) : (method == :average ? 0. : 2.)) #0. is for the average update, 2. is when reducing the τ matrix (:contraction)
    for k in neighbors(g, c) 
        if k ≠ v
            D[v,k] = (distance(D,i,k) + distance(D,k,j) - reduction)/2
            D[k,v] = D[v,k]
        end
    end
    D[v,v] = 0.
end

@inline function distance(D::AbstractMatrix, i, j)
    D[max(i,j), min(i,j)] # nodes always have up to date average distance wrt nodes with a smaller index
end
@inline function distance(D::AbstractMatrix, vect)
    i,j = vect
    distance(D, i, j) 
end
