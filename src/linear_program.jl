using JuMP, Combinatorics
export lp_relaxation_x_tau, weak_buneman_oracle

function lp_exact_x_tau(g, D::Matrix, c; relax = true) # D is a distance matrix (2n-2 x 2n-2), g is a graph containing a star with n branches and internal node c = n+1 
    n = degree(g,c)
    TAXA = sort(collect(neighbors(g,c)))
    LENGTHS = 2:n-1
    model = Model(OPTIMIZER, add_bridges = false)
    set_string_names_on_creation(model, false)
    @variable(model, _x[i=TAXA, j=TAXA, l=LENGTHS; i<j], Bin)
    @expression(model, x[i=TAXA, j=TAXA, l=LENGTHS; i≠j], _x[minmax(i,j)...,l])
    @expression(model, τ[i=TAXA, j=TAXA; i≠j], sum(l*_x[minmax(i,j)...,l] for l in LENGTHS))
    @variable(model, 1<=s[i=TAXA[1:end-1], j=TAXA, k=TAXA; i≠j && j≠k]<=n-2, Int)
    bun_idxs = vcat([(j,p,q) for (j,p,q) in combinations(TAXA, 3) if j > TAXA[1]], [(p,j,q) for (j,p,q) in combinations(TAXA, 3) if j > TAXA[1]], [(q,j,p) for (j,p,q) in combinations(TAXA, 3) if j > TAXA[1]])
    @variable(model, y[bun_idxs], Bin)

    @constraints model begin
        c_one_length[i=TAXA, j=TAXA; i<j], 
            sum(x[i,j,l] for l in LENGTHS) == 1
        c_kraft[i=TAXA], 
            sum(sum(exp2(-l)*x[i,j,l] for l in LENGTHS) for j in TAXA if j ≠ i)  == 0.5
        c_manifold, 
            sum(l*exp2(-l)*x[i,j,l] for i in TAXA, j in TAXA, l in LENGTHS if i≠j) == (2n-3)
        c_str_triangular[i=TAXA[1:end-1], j=TAXA, k=TAXA; i≠j && j≠k && k > i],
            τ[j,k] + τ[i,j] == 2*s[i,j,k] + τ[i,k]
        c_bun_sum[j in TAXA[2:end-2], p in TAXA[1:end-1], q in TAXA; p>j && q>p],
            y[(j,p,q)] + y[(p,j,q)] + y[(q,j,p)] == 1            

    end
    for (j, p, q) in combinations(TAXA, 3)
        if j == first(TAXA)
            continue
        end
        @constraint(model, s[first(TAXA),p,j] >= 1 + y[(j,minmax(p,q)...)])
        @constraint(model, s[first(TAXA),p,q] >= 1 + y[(q,minmax(p,j)...)])
        @constraint(model, s[first(TAXA),j,p] >= 1 + y[(p,minmax(j,q)...)])
        @constraint(model, s[first(TAXA),j,q] >= 1 + y[(q,minmax(j,p)...)])
        @constraint(model, s[first(TAXA),q,j] >= 1 + y[(j,minmax(p,q)...)])
        @constraint(model, s[first(TAXA),q,p] >= 1 + y[(p,minmax(q,j)...)])
    end
    @objective(model, Min, sum(D[i,j]*sum(exp2(-l)*x[i,j,l] for l in LENGTHS) for i in TAXA, j in TAXA if i≠j))
    if relax
        relax_integrality(model)
    end
    return model
end


function lp_relaxation_x_tau(g, D::Matrix, c; relax = true) # D is a distance matrix (2n-2 x 2n-2), g is a graph containing a star with n branches and internal node c = n+1 
    n = degree(g,c)
    TAXA = sort(collect(neighbors(g,c)))
    LENGTHS = 2:n-1
    model = Model(OPTIMIZER, add_bridges = false)
    set_string_names_on_creation(model, false)
    @variable(model, _x[i=TAXA, j=TAXA, l=LENGTHS; i<j], Bin)
    @expression(model, x[i=TAXA, j=TAXA, l=LENGTHS; i≠j], _x[minmax(i,j)...,l])
    @expression(model, τ[i=TAXA, j=TAXA; i≠j], sum(l*_x[minmax(i,j)...,l] for l in LENGTHS))
    @constraints model begin
        c_one_length[i=TAXA, j=TAXA; i<j], 
            sum(x[i,j,l] for l in LENGTHS) == 1
        c_kraft[i=TAXA], 
            sum(sum(exp2(-l)*x[i,j,l] for l in LENGTHS) for j in TAXA if j ≠ i)  == 0.5
        c_manifold, 
            sum(l*exp2(-l)*x[i,j,l] for i in TAXA, j in TAXA, l in LENGTHS if i≠j) == (2n-3)
    end
    @objective(model, Min, sum(D[i,j]*sum(exp2(-l)*x[i,j,l] for l in LENGTHS) for i in TAXA, j in TAXA if i≠j))
    if relax
        relax_integrality(model)
    end
    return model
end
