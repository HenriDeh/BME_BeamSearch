using GLMakie, GraphMakie, Colors, NetworkLayout
export plot_ubt
GLMakie.activate!(inline=false)

function plot_ubt(g; center = nothing, subleaves = [], kwargs...)
    nleaves = (nv(g)+2)÷2
    labels = string.(1:nv(g)) #[string.(1:nleaves); fill("", nleaves-2)]
    if center isa Int
        subleavesset = Set(subleaves)
        if !isempty(subleavesset)
            new_inners = Set(neighbors(g,center))
            inners = Set{Int}()
            while !isempty(new_inners)
                n = pop!(new_inners)
                push!(inners, n)
                union!(new_inners, setdiff(neighbors(g,n), subleavesset, inners))
            end
        end
        colors = [n in subleavesset ? colorant"lightgreen" : n == center ? colorant"yellow" : n in inners ? colorant"orange" : colorant"turquoise" for n in 1:nv(g)] 
    else
        colors = [fill(colorant"lightgreen", nleaves); fill(colorant"lightgrey", nleaves-2)] 
    end
    f = Figure(size = (1000,1000))
    ax = Axis(f[1,1], title = get(kwargs, :title, ""))
    p = graphplot!(ax, g, ilabels = labels, node_size = fill(20, nv(g)), node_color = colors, kwargs...)
    deregister_interaction!(ax, :rectanglezoom)
    register_interaction!(ax, :nodedrag, NodeDrag(p))
    register_interaction!(ax, :nodehover, NodeHoverHighlight(p))
    hidedecorations!(ax); 
    hidespines!(ax);
    return f
end


