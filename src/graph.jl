using Reexport
using Combinatorics, InvertedIndices, OffsetArrays, StatsBase
import DataStructures.Stack
@reexport using Graphs

export UBT_from_PLM, random_subtree, create_random_UBT, path_length_matrix, collaspe_to_inner_star, star_graph, tree_length

function tree_length(plm::Matrix, D)
    sum(D[i,j]*2.0^-plm[i,j] for (i,j) in combinations(1:size(plm,1),2), init = 0.)*2
end

function find_cherry(D::Matrix, neighbors)
    for i in neighbors
        for j in neighbors
            i < j || continue
            D[i,j] == 2 && return (i,j)
        end
    end
    return nothing
end

star_graph(D::Matrix) = star_graph(size(D,1))

function star_graph(n::Int)
    g = Graph(n+1)
    c = n+1
    for i in 1:c-1
        add_edge!(g, i, c)
    end
    return g, c
end

"""
    UBT_from_PLM(PLM::Matrix)

Instantiate the graph object of the UBT represented by the path-length matrix PLM.

"""
function UBT_from_PLM(PLM::Matrix)
    n = size(PLM,1)
    g = Graph(n+1)
    c = n+1
    for i in 1:n 
        add_edge!(g, i, c)
    end
    _PLM = extend_distance(PLM)
    while degree(g, c) > 3
        (i,j) = find_cherry(_PLM, neighbors(g,c))
        add_vertex!(g)
        v = nv(g)
        rem_edge!(g, c, i)
        rem_edge!(g, c, j)
        add_edge!(g, v, i)
        add_edge!(g, v, j)
        add_edge!(g, v, c)
        for k in neighbors(g,c)
            if k == v
                _PLM[v,v] = 0.
            else
                if (i,k,j) == (3,2,4)
                    @assert _PLM[i,k] == _PLM[j, k] == _PLM[k,i] == _PLM[k,j]
                end
                _PLM[v,k] = _PLM[i,k] - 1
                _PLM[k,v] = _PLM[v,k]
            end
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

function collaspe_to_inner_star(_g, subtree, D::Matrix) #D is the distance matrix of the leaves, K is number of star branches
    g = copy(_g)
    selected = subtree.inner_nodes
    leaves = subtree.taxa_leaves
    candidates = subtree.prune_leaves
    ## make star
    add_vertex!(g) #center of the star
    c = nv(g)
    for i in union(leaves, candidates)
        add_edge!(g,c,i) #create the star branches
    end
    node_map = rem_vertices!(g, collect(selected)) #remove inner nodes of subgraph
    c = findfirst(==(c), node_map)
    D2 = extend_distance(D)
    parents = Stack{Pair{Tuple{Int,Int},Int}}()
    visited = Set{Int}(vcat([c], neighbors(g,c)))
    remaining = Set{Int}() # set of nodes that have not yet been contracted. 
    for v in BFSIterator(g, c) # recursively obtain childs-parent, starting from c to search outwards towards the leaves of the UBT
        v == c && continue
        if degree(g, v) == 1
            push!(remaining, v)
            continue
        end
        c1, c2 = [n for n in neighbors(g,v) if n ∉ visited] #avoid backtracking
        push!(visited, c1)
        push!(visited, c2)
        push!(parents, (c1, c2) => v)
    end
    while !isempty(parents)
        (c1, c2), p = pop!(parents)
        pop!(remaining, c1)
        pop!(remaining, c2)
        for v in remaining
            D2[v, p] = D2[p, v] = (D2[c1, v] + D2[c2, v])/2# - D2[c1,c2])/2 #contraction of c1 and c2 into their parent.
        end
        push!(remaining, p)
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

#leavesX contains the depth of each leave of the subtree X
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

function is_ubt(g)
    c = counter(degree(g))
    return c[1] == length(leaves(g)) && c[3] == length(inner_nodes(g))
end

#Catanzaro, D., & Pesenti, R. (2019). Enumerating vertices of the balanced minimum evolution polytope. Computers & Operations Research, 109, 209‑217. https://doi.org/10.1016/j.cor.2019.05.001
function path_length_matrix(ubt::Graph)
    @assert is_ubt(ubt) "Graph is not an UBT"
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

function compute_splits(g::Graph) 
    if length(vertices(g)) == length(edges(g)) + 1
        compute_splits(g, length(leaves(g))+1)
    else
        root = sum(degree(g) .== 1) + 1
        compute_splits(g, root)
    end
end
# root must be |Γ|+1, usually the center of the star_graph
function compute_splits(g::Graph, root)
    taxa = Set(collect(1:root-1))

    subtrees = Dict{Int, Set{Int}}()
    splits = Set{Tuple{Set{Int},Set{Int}}}()
    stack = Stack{Tuple{Int,Int}}()
    push!(stack, (0, root)) # root of postorder traversal, 0 is dummy parent
    visited = Set{Int}() # to not revisit parents as g is undirected

    # Stacks for processing internal nodes after children are visited
    postorder_stack = Stack{Tuple{Int,Int}}()
    while !isempty(stack)
        p, n = pop!(stack)
        push!(postorder_stack, (p,n))
        push!(visited, n)

        for c in neighbors(g, n)
            c in visited && continue
            push!(stack, (n,c))
        end
    end
    # Process nodes in postorder and compute the subtrees of each node
    while !isempty(postorder_stack)
        (p, n) = pop!(postorder_stack)
        if degree(g, n) == 1
            subtrees[n] = Set([n])
        else
            subtrees[n] = union([subtrees[c] for c in neighbors(g, n) if c != p]...)
        end
        p == 0 && continue
        a,b = subtrees[n], setdiff(taxa, subtrees[n])
        push!(splits, (a,b))
        push!(splits, (b,a))
    end
    return splits
end

are_same_topologies(g1::Graph, g2::Graph) = are_same_topologies(compute_splits(g1), compute_splits(g2))

function are_same_topologies(splitsT1, splitsT2)
    if length(splitsT1) != length(splitsT2) 
        @error "Trees do not have same number of splits"
    end
    return all(in(splitsT2), splitsT1)
end