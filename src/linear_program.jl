using JuMP, Combinatorics
export MIP_reduced, weak_buneman_oracle

function MIP_complete(g, D::Matrix, c; relax = true, semi_relax = false) # D is a distance matrix (2n-2 x 2n-2), g is a graph containing a star with n branches and internal node c = n+1 
    n = degree(g,c)
    TAXA = sort(collect(neighbors(g,c)))
    taxon1 = first(TAXA)
    TAXA2n = TAXA[2:end]
    LENGTHS = 2:n-1
    model = Model(OPTIMIZER, add_bridges = false)
    set_string_names_on_creation(model, false)
    if semi_relax
        @variable(model, 0<=_x[i=TAXA, j=TAXA, l=LENGTHS; i<j]<=1)
    else
        @variable(model, _x[i=TAXA, j=TAXA, l=LENGTHS; i<j], Bin)
    end
    @variable(model, _y[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p < q], Bin)
    
    @expression(model, y[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p ≠ q], _y[j, minmax(p,q)...])
    @expression(model, x[i=TAXA, j=TAXA, l=LENGTHS; i≠j], _x[minmax(i,j)...,l])
    @expression(model, τ[i=TAXA, j=TAXA; i≠j], sum(l*_x[minmax(i,j)...,l] for l in LENGTHS))

    @constraints model begin
        c_one_length[i=TAXA, j=TAXA; i<j], 
            sum(x[i,j,l] for l in LENGTHS) == 1
        c_kraft[i=TAXA], 
            sum(sum(exp2(-l)*x[i,j,l] for l in LENGTHS) for j in TAXA if j ≠ i)  == 0.5
        c_manifold, 
            sum(l*exp2(-l)*x[i,j,l] for i in TAXA, j in TAXA, l in LENGTHS if i≠j) == (2n-3)
        
        c_str_triangular[i=TAXA2n, j=TAXA2n; i ≠ j],
            τ[taxon1, i] + τ[taxon1, j] - τ[i,j] >= 2
        c_str_triangular2[i=TAXA2n, j=TAXA2n; i ≠ j],
            τ[taxon1, i] + τ[i,j] - τ[taxon1,j] >= 2
        
        c_bun_sum[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p ≠ q],
            y[j,p,q] + y[p,j,q] + y[q,j,p] == 1            
        c_buneman[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p ≠ q],
            τ[taxon1,j] + τ[p,q] >= 2(1 - y[q,j,p]) + τ[taxon1,p] + τ[j,q] - (2n - 2)y[j,p,q]
    end

    @objective(model, Min, sum(D[i,j]*sum(exp2(-l)*x[i,j,l] for l in LENGTHS) for i in TAXA, j in TAXA if i≠j))
    if relax
        relax_integrality(model)
    end
    return model
end


function MIP_reduced(g, D::Matrix, c; relax = true) # D is a distance matrix (2n-2 x 2n-2), g is a graph containing a star with n branches and internal node c = n+1 
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
