using JuMP, Combinatorics
export BMEP_MILP

function BMEP_MILP(D::Matrix; relax = true, semi_relax = false, complete = true, init_tau = fill(2, size(D)), scale = false, max_length = true) # D is a distance matrix (n x n)
    n = size(D,1)
    TAXA = 1:n
    taxon1 = first(TAXA)
    TAXA2n = TAXA[2:end]
    LENGTHS = 2:(max_length ? min(n-1, ceil(log2(n-1)^2)) : n-1)
    model = Model(OPTIMIZER, add_bridges = false)
    set_string_names_on_creation(model, false)
    if semi_relax || relax #option to only optimize on binary y's and let x be fractional
        @variable(model, 0 <= _x[i=TAXA, j=TAXA, l=LENGTHS; i<j] <= 1, start = init_tau[i,j] == l ? 1 : 0)
    else
        @variable(model, _x[i=TAXA, j=TAXA, l=LENGTHS; i<j], Bin, start = init_tau[i,j] == l ? 1 : 0)
    end
    if complete
        if relax 
            @variable(model, 0 <= _y[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p < q] <= 1)
        else
            @variable(model, _y[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p < q], Bin)
        end
        @expression(model, y[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p ≠ q], _y[j, minmax(p,q)...])
    end
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
        if complete
            @constraints model begin
                c_str_triangular[i=TAXA2n, j=TAXA2n; i ≠ j],
                    τ[taxon1, i] + τ[taxon1, j] - τ[i,j] >= 2
                c_str_triangular2[i=TAXA2n, j=TAXA2n; i ≠ j],
                    τ[taxon1, i] + τ[i,j] - τ[taxon1,j] >= 2
                c_bun_sum[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p ≠ q],
                    y[j,p,q] + y[p,j,q] + y[q,j,p] == 1            
                c_buneman[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p ≠ q],
                    τ[taxon1,j] + τ[p,q] >= 2(1 - y[q,j,p]) + τ[taxon1,p] + τ[j,q] - (2n - 2)y[j,p,q]
        end
    end
    if scale
        @objective(model, Min, sum(D[i,j]*sum(exp2(n-1-(l-1))*x[i,j,l] for l in LENGTHS) for i in TAXA, j in TAXA if i<j))
    else
        @objective(model, Min, sum(D[i,j]*sum(exp2(-(l-1))*x[i,j,l] for l in LENGTHS) for i in TAXA, j in TAXA if i<j))
    end
    return model
end

function set_distances!(model, D; scale = false, max_length = true)
    n = size(D,1)
    TAXA = 1:n
    LENGTHS = 2:(max_length ? min(n-1, ceil(log2(n-1)^2)) : n-1)
    x = model[:x]
    if scale
        @objective(model, Min, sum(D[i,j]*sum(exp2(n-1-(l-1))*x[i,j,l] for l in LENGTHS) for i in TAXA, j in TAXA if i<j))
    else
        @objective(model, Min, sum(D[i,j]*sum(exp2(-(l-1))*x[i,j,l] for l in LENGTHS) for i in TAXA, j in TAXA if i<j))
    end
end
