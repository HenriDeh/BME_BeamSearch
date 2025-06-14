using FastME_jll, DelimitedFiles
using Phylo: parsenewick, traversal, isleaf, isroot, getparent, preorder
import DataStructures.Stack
export fastme_local_search, read_distance_matrix, ubtgraph_from_nwk

"""
    fastme_help()

Show help for the FastME command-line tool.
"""
fastme_help() = run(`$(fastme()) -h`)

"""
    ubtgraph_from_nwk(path::String)

Parse a Newick file at `path` and return the corresponding UBT graph.
"""
function ubtgraph_from_nwk(path::String)
    t = open(path) do f
        ln = readline(f)
        l = replace(ln, r":-?\d+\.\d+," => ",", r":-?\d+\.\d+\)" => ")")
        #l = replace(l, r"\) [^:]*\[[^:]*:" => ") :") # remove branch lengths, remove buggy encoding node names
        parsenewick(l)
    end
    ubtgraph_from_nwk(t)
end

"""
    ubtgraph_from_nwk(t)

Convert a parsed Newick tree `t` to a UBT graph.
"""
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

"""
    fastme_local_search(D::Matrix, g::Graph)

Run FastME local search starting from graph `g` and distance matrix `D`.

# Arguments
- `D::Matrix`: Distance matrix.
- `g::Graph`: Initial tree (UBT graph) to start the search from.
"""
function fastme_local_search(D::Matrix, g::Graph)
    tmp_tree = "testtree.nwk"
    open(tmp_tree, "w") do f
        write(f, ubt_to_nwk(g))
    end

    gspr = try 
        fastme_local_search(D, tmp_tree, inittree = true)
    finally
        rm(tmp_tree, force = true)
    end
    return gspr
end

"""
    fastme_local_search(D::Matrix, path::String, tree_path::String; inittree = false, verbose = true, spr = true, methods = ["b","o","i","n","u"])

Run FastME local search with distance matrix `D` and initial tree at `tree_path`.

# Arguments
- `D::Matrix`: Distance matrix.
- `path::String`: Path to write the distance matrix.
- `tree_path::String`: Path to write/read the initial tree.
- `inittree`: Use initial tree (`tree_path`) if true.
- `verbose`: Print progress if true.
- `spr`: Use subtree pruning and regrafting if true.
- `methods`: List of FastME initialization methods to try.
"""
function fastme_local_search(D::Matrix, path::String, tree_path::String; inittree = false, verbose = true, spr = true, methods = ["b","o","i","n","u"])
    D_to_txt(path, D)
    fastme_local_search(path, tree_path; methods, inittree, verbose, spr)
end

"""
    fastme_local_search(D::Matrix, tree_path::String; inittree = false, verbose = true, spr = true, methods = ["b","o","i","n","u"])

Run FastME local search with distance matrix `D` and output to `tree_path`.

# Arguments
- `D::Matrix`: Distance matrix.
- `tree_path::String`: Path to output tree file.
- `inittree`: Use initial tree if true.
- `verbose`: Print progress if true.
- `spr`: Use subtree pruning and regrafting if true.
- `methods`: List of FastME initialization methods to try.
"""
function fastme_local_search(D::Matrix, tree_path::String; inittree = false, verbose = true, spr = true, methods = ["b","o","i","n","u"])
    D_to_txt("tmpDfm.txt", D)
    try 
        return fastme_local_search("tmpDfm.txt", tree_path; methods, inittree, verbose, spr)
    finally
        rm("tmpDfm.txt", force = true)
    end
end

"""
    fastme_local_search(D::Matrix, τ::Matrix, path::String, tree_path::String)

Run FastME local search with distance matrix `D`, path-length matrix `τ`, and output to `tree_path`.
"""
function fastme_local_search(D::Matrix, τ::Matrix, path::String, tree_path::String)
    D_to_txt(path,D)
    nwk = PLM_to_nwk(τ)
    open(tree_path, "w") do f
        println(f, nwk)
    end
    fastme_local_search(path, tree_path, inittree = true)
end

"""
    fastme_local_search(path::String, tree_path::String; inittree = false, methods = ["b","o","i","n","u"], verbose = true, spr = true)

Run FastME local search with distance matrix at `path` and output to `tree_path`.

# Arguments
- `path::String`: Path to distance matrix file.
- `tree_path::String`: Path to output tree file.
- `inittree`: Use initial tree if true.
- `methods`: List of FastME initialization methods to try.
- `verbose`: Print progress if true.
- `spr`: Use subtree pruning and regrafting if true.
"""
function fastme_local_search(path::String, tree_path::String; inittree = false, methods = ["b","o","i","n","u"], verbose = true, spr = true)
    if inittree
        run(pipeline(`$(fastmeMP()) --spr -i $path -T $(Threads.nthreads()) -o tmp.nwk -w BalLS -m u -u $tree_path`, stdout=devnull))
        mv("tmp.nwk", tree_path, force = true)
    else
        @assert all(in(["b","o","i","n","u"]), methods)
        D = read_distance_matrix(path)
        best = Inf
        for method in methods
            if spr
                run(pipeline(`$(fastmeMP()) -i $path --spr -T $(Threads.nthreads()) -o tmp.nwk -w BalLS -m $method`, stdout=devnull))
            else
                run(pipeline(`$(fastmeMP()) -i $path -T $(Threads.nthreads()) -o tmp.nwk -w BalLS -m $method`, stdout=devnull))
            end
            try
                g = ubtgraph_from_nwk("tmp.nwk")
                tl = tree_length(path_length_matrix(g), D)
                if verbose
                    println("FastME$(spr ? "+spr" : "") init method $method yielded tree length $tl")
                end
                if tl < best
                    best = tl
                    mv("tmp.nwk", tree_path, force = true)
                else
                    rm("tmp.nwk")
                end
            catch e #some -b methods of FastME may yield incorrect newicks on Windows
                if verbose
                    println("Method $method errored")
                end
                rm("tmp.nwk")
                #throw(e)
                continue
            end
        end

    end
    return ubtgraph_from_nwk(tree_path)
end

"""
    read_distance_matrix(path::String)

Read a distance matrix from a file at `path`.
"""
function read_distance_matrix(path::String)
    D = open(path) do file
        M = readdlm(file, header = true)[1]
        return M[:,2:end]
    end 
end

"""
    ubt_to_nwk(g::Graph)

Convert a UBT graph `g` to a Newick string.
"""
function ubt_to_nwk(g::Graph)
    n = length(leaves(g))
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
#        nwk *= "$(reverse(string(current)))"
        w = reverse(string(1.0)) # we give each edge a weight of 1 otherwise FastME thinks the tree has a length of 0.
        if current <= n
            nwk *= "$w:$(reverse(string(current)))"
        else
            nwk *= "$w:"
        end


        for n in outneighbors(bf, current)
            push!(S, n)
            push!(depth_stack, depth+1)
        end
        prev_depth = depth
    end
    nwk *= "("^(prev_depth-pop!(depth_stack))
    return reverse(nwk)
end

"""
    PLM_to_nwk(τ::Matrix)

Convert a path-length matrix `τ` to a Newick string.
"""
function PLM_to_nwk(τ::Matrix)
    nwk = ubt_to_nwk(UBT_from_PLM(τ))
end
