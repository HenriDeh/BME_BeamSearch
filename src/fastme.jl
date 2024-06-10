using FastME_jll, DelimitedFiles
using Phylo: parsenewick, traversal, isleaf, isroot, getparent, preorder
import DataStructures.Stack
export fastme_local_search, read_distance_matrix, ubtgraph_from_nwk

fastme_help() = run(`$(fastme()) -h`)

function ubtgraph_from_nwk(path::String)
    t = open(path) do f
        l = replace(readline(f), r":-?\d+\.\d+," => ",", r":-?\d+\.\d+\)" => ")")
        parsenewick(l)
    end
    ubtgraph_from_nwk(t)
end

function ubtgraph_from_nwk(t)
    trav = traversal(t,preorder)
    g = Graph()
    add_vertices!(g, (length(t.nodes)+2)÷2)
    nodemap = Dict{String, Int}()
    for node in trav
        if !isleaf(t, node)
            add_vertex!(g)
            nodemap[node.name] = nv(g)
        else
            nodemap[node.name] = parse(Int, node.name)
        end
        if !isroot(t, node)
            p = nodemap[getparent(t,node).name]
            add_edge!(g, p, nodemap[node.name])

        end
    end
    return g
end

function fastme_local_search(path::String, tree_path::String; inittree = false)
    if inittree
        run(pipeline(`$(fastmeMP()) -i $path --spr -T $(Threads.nthreads()) -o tmp.nwk -w none -u $tree_path`, stdout=devnull))
    else
        run(pipeline(`$(fastmeMP()) -i $path --spr -T $(Threads.nthreads()) -o tmp.nwk -w none -m B`, stdout=devnull))
    end
    mv("tmp.nwk", tree_path, force = true)
    g = ubtgraph_from_nwk(tree_path)
    return g
end

function read_distance_matrix(path::String)
    D = open(path) do file
        M = readdlm(file, header = true)[1]
        return M[:,2:end]
    end 
end

function ubt_to_nwk(g::Graph)
    nwk = ";"
    current = (nv(g)+2)÷2+1
    bf = bfs_tree(g,current)
    S = Stack{Int}()
    push!(S, current)
    depth_stack = Stack{Int}()
    push!(depth_stack, 0)
    push!(depth_stack, 0)
    prev_depth = 0
    while !isempty(S)

        current = pop!(S)
        depth = pop!(depth_stack)

        if prev_depth < depth
            nwk *= ")"
        elseif prev_depth > depth
            nwk *= "("^(prev_depth-depth)
            nwk *= ","
        elseif length(nwk) > 1
            nwk *= ","
        end
        nwk *= "$(reverse(string(current)))"

        for n in outneighbors(bf, current)
            push!(S, n)
            push!(depth_stack, depth+1)
        end
        prev_depth = depth
    end
    nwk *= "("^(prev_depth-pop!(depth_stack))
    return reverse(nwk)
end
