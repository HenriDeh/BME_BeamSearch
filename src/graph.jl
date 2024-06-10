using Reexport
using Combinatorics, InvertedIndices, OffsetArrays, StatsBase
import DataStructures.Stack
@reexport using Graphs

export make_graph, random_subtree, create_random_UBT, path_length_matrix, collaspe_to_inner_star

function find_cherry(D::Dict)
    for i in keys(D)
        for j in keys(D[i])
            D[i][j] == 2 && return (i,j)
        end
    end
    return nothing
end

"""
    make_graph(Tau::Matrix)

Instantiate the graph object of the UBT represented by the path-length matrix Tau.

"""
function make_graph(Tau::Matrix)
    n = size(Tau,1)
    g = Graph(n+1)
    for i in 1:n 
        add_edge!(g, i, n+1)
    end
    distances = Dict(i => Dict(j => Tau[i,j] for j in 1:n) for i in 1:n)
    for v in n+2:(2n-2)
        (i,j) = find_cherry(distances)
        add_vertex!(g)
        rem_edge!(g, n+1, i)
        rem_edge!(g, n+1, j)
        add_edge!(g, v, i)
        add_edge!(g, v, j)
        add_edge!(g, v, n+1)
        distances[v] = Dict(k => d-1 for (k,d) in distances[i])
        pop!(distances, i)
        pop!(distances, j)
        for k in keys(distances)
            distances[k][v] = distances[k][i] - 1 
            pop!(distances[k], i)
            pop!(distances[k], j)
        end
    end
    return g
end

leaves(g) = 1:(nv(g)÷2+1)
inner_nodes(g) = (nv(g)÷2+2):nv(g)

"""
    random_subtree(g, K [, D])

Samples a subtree of g that is also a UBT with K leaves. 
Returns `sg, node_map` where `sg` is the subgraph object and `node_map` maps the numbering of the 
subgraph nodes to the numbering of the original graph `g`.

If D, a symmetric distance matrix between the leaves of `g`, is provided, returns 
`sg, node_map, D2, leave_map` where D2 is the distance matrix between the leaves of 
`sg` computed with the recursive average distance. `leave_map` maps the indices of 
that matrix to the original nodes of `g`. 
"""
function random_subtree(g, K)
    start = rand(inner_nodes(g))
    selected = Set([start])
    leaves = Set(filter(n->degree(g,n) == 1, neighbors(g,start)))
    candidates = Set{Int}(filter(n->degree(g,n) == 3, neighbors(g,start)))
    println(selected)
    println(candidates)
    while length(selected) < K-2
        println()
        new_node = rand(candidates) 
        println(new_node)
        push!(selected, new_node)
        pop!(candidates, new_node) #remove it from candidates
        union!(candidates, setdiff(filter(n->degree(g,n) == 3, neighbors(g,new_node)), selected)) #add its neighbors in candidates
        union!(leaves, filter(n->degree(g,n) == 1, neighbors(g,new_node)))
        println(candidates)
        println(leaves)
    end
    println(selected)
    return induced_subgraph(g, collect(union(selected, leaves, candidates)))
end

function collaspe_to_inner_star(_g, K, D::Matrix; neighborhood) #D is the distance matrix of the leaves, K is number of star branches
    g = copy(_g)
    selected, leaves, candidates = neighborhood(g, K)
    nodes = union(selected, leaves, candidates)
    subtree_depths = Dict{Int, Dict{Int,Int}}() #std[i] is a dictionary with all nodes of the substree rooted in i and their associated depths in this subtree.
    for ca in candidates
        subtree_depths[ca] = find_subtree_depths(ca, g, nodes)
    end
    for l in leaves
        subtree_depths[l] = Dict(l => 0)
    end
    ## make star
    add_vertex!(g) #center of the star
    c = nv(g)
    for i in union(leaves, candidates)
        add_edge!(g,c,i) #create the star branches
    end
    node_map = rem_vertices!(g, collect(selected)) #remove inner nodes of subgraph
    c = findfirst(==(c), node_map)
    update_weights!(neighborhood, node_map)
    D2 = extend_distance(D)
    for (i,j) in combinations(neighbors(g,c),2)
        D2[i,j] = average_dist(D, subtree_depths[node_map[i]], subtree_depths[node_map[j]])
        D2[j,i] = D2[i,j]
    end
    return g, D2, c, node_map
end

function find_subtree_depths(c, g, subtree_nodes)
    Γ = 1:(nv(g)+2)÷2
    to_expand = Stack{Tuple{Int,Int}}() #(node, depth)
    push!(to_expand, (c, 0))
    leaves_depths = Dict{Int,Int}()
    while !isempty(to_expand)
        k, depth = pop!(to_expand)
        push!(subtree_nodes, k)
        l, m = setdiff(neighbors(g, k), subtree_nodes)
        if l ∉ Γ
            push!(to_expand, (l, depth + 1))
        else
            push!(leaves_depths, l => depth + 1)
        end
        if m ∉ Γ
            push!(to_expand, (m, depth + 1))
        else
            push!(leaves_depths, m => depth + 1)
        end        
    end
    return leaves_depths
end

function average_dist(D::Matrix, leavesA, leavesB)
    @inbounds sum(D[l,m]*2.0^-(leavesA[l]+leavesB[m]) for (l,m) in Iterators.product(keys(leavesA), keys(leavesB)))
end

function create_random_UBT(n)
    g = Graph(n+1)
    for i in 1:n 
        add_edge!(g, i, n+1)
    end
    candidates = Set(collect(1:n))
    while nv(g) < 2n -2
        i = rand(candidates)
        pop!(candidates, i)
        j = rand(candidates)
        pop!(candidates, j)
        rem_edge!(g,n+1, i)
        rem_edge!(g,n+1, j)
        add_vertex!(g)
        add_edge!(g,i, nv(g))
        add_edge!(g,j, nv(g))
        add_edge!(g,nv(g),n+1)
        push!(candidates, nv(g))
    end
    return g
end

#Catanzaro, D., & Pesenti, R. (2019). Enumerating vertices of the balanced minimum evolution polytope. Computers & Operations Research, 109, 209‑217. https://doi.org/10.1016/j.cor.2019.05.001
function path_length_matrix(ubt::Graph)
    s = prufer_encode(ubt)
    ir = inner_nodes(ubt)
    k = length(ir)
    M = OffsetMatrix(zeros(Int, k, 6), ir, 1:6)
    τ = zeros(Int, k+2, k+2)
    for i in ir
        M[i,1:3] .= neighbors(ubt, i)
    end
    for i in leaves(ubt)
        M[s[i], 4:6] .= [1, i, 0]
        current_node = s[i]
        while true
            if M[current_node, 6] < 3
                M[current_node, 6] += 1
                visited_node = M[current_node, M[current_node, 6]]
                if visited_node != M[current_node, 5]
                    if visited_node in ir
                        M[visited_node, 4] = M[current_node, 4] + 1
                        M[visited_node, 5] = current_node
                        M[visited_node, 6] = 0
                        current_node = visited_node
                    else
                       τ[i, visited_node] = M[current_node, 4] + 1
                       τ[visited_node, i] = τ[i, visited_node]
                    end
                end
            else
                if current_node == s[i]
                    break
                else
                    current_node = M[current_node, 5]
                end
            end
        end
    end
    return τ
end
