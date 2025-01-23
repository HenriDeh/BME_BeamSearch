using JuMP, Combinatorics
export BMEP_MILP

function BMEP_MILP(D::Matrix; relax = true, triangular = true, buneman = true, init_tau = fill(2, size(D)), scale = false, max_length = false) # D is a distance matrix (n x n)
    n = size(D,1)
    TAXA = 1:n
    taxon1 = first(TAXA)
    TAXA2n = TAXA[2:end]
    Lmax = max_length ? min(n-1, ceil(log2(n-1)^2)) : n-1
    LENGTHS = 2:Lmax
    model = Model(OPTIMIZER, add_bridges = false)
    set_string_names_on_creation(model, false)
    @expression model n n
    if relax
        @variable(model, 0 <= _x[i=TAXA, j=TAXA, l=LENGTHS; i<j] <= 1, start = init_tau[i,j] == l ? 1 : 0)
    else
        @variable(model, _x[i=TAXA, j=TAXA, l=LENGTHS; i<j], Bin, start = init_tau[i,j] == l ? 1 : 0)
    end
    if buneman
        if relax 
            @variable(model, 0 <= _y[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p < q] <= 1)
        else
            @variable(model, _y[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p < q], Bin)
        end
        @expression(model, y[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p ≠ q], _y[j, minmax(p,q)...])
    end
    @expression(model, x[i=TAXA, j=TAXA, l=LENGTHS; i≠j], _x[minmax(i,j)...,l])
    @expression(model, τ[i=TAXA, j=TAXA; i≠j], sum(l*_x[minmax(i,j)...,l] for l in LENGTHS))
    # if n > 3
    #     @constraint model c_cherry[i = TAXA] sum(x[i,j,2] for j in TAXA if j ≠ i) <= 1
    # end
    @constraints model begin
        c_one_length[i=TAXA, j=TAXA; i<j], 
            sum(x[i,j,l] for l in LENGTHS) == 1
        c_kraft[i=TAXA], 
            sum(sum(exp2(-l)*x[i,j,l] for l in LENGTHS) for j in TAXA if j ≠ i)  == 0.5
        c_manifold, 
            sum(l*exp2(-l)*x[i,j,l] for i in TAXA, j in TAXA, l in LENGTHS if i≠j) == (2n-3)
        end
        if triangular
            @constraints model begin
                c_str_triangular[i=TAXA2n, j=TAXA2n; i ≠ j],
                    τ[taxon1, i] + τ[taxon1, j] - τ[i,j] >= 2
                c_str_triangular2[i=TAXA2n, j=TAXA2n; i ≠ j],
                    τ[taxon1, i] + τ[i,j] - τ[taxon1,j] >= 2
            end
        end
        if buneman
            @constraints model begin
                c_bun_sum[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p ≠ q],
                    y[j,p,q] + y[p,j,q] + y[q,j,p] == 1            
                c_buneman[j = TAXA2n, p = TAXA2n, q = TAXA2n; j ≠ p && j ≠ q && p ≠ q],
                    τ[taxon1,j] + τ[p,q] >= 2(1 - y[q,j,p]) + τ[taxon1,p] + τ[j,q] - (2n - 2)y[j,p,q]
        end
    end
    @expression(model, tree_length, sum(D[i,j]*sum(exp2(-(l-1))*x[i,j,l] for l in LENGTHS) for i in TAXA, j in TAXA if i<j))
    if scale
        @objective(model, Min, sum(D[i,j]*sum(exp2(Lmax-(l-1))*x[i,j,l] for l in LENGTHS) for i in TAXA, j in TAXA if i<j))
    else
        @objective(model, Min, tree_length)
    end
    return model
end

function set_distances!(model, D; scale = false, max_length = true)
    n = size(D,1)
    TAXA = 1:n
    Lmax = max_length ? min(n-1, ceil(log2(n-1)^2)) : n-1
    LENGTHS = 2:Lmax
    x = model[:x]
    unregister(model, :tree_length)
    @expression(model, tree_length, sum(D[i,j]*sum(exp2(-(l-1))*x[i,j,l] for l in LENGTHS) for i in TAXA, j in TAXA if i<j))
    if scale
        @objective(model, Min, sum(D[i,j]*sum(exp2(Lmax-(l-1))*x[i,j,l] for l in LENGTHS) for i in TAXA, j in TAXA if i<j))
    else
        @objective(model, Min, tree_length)
    end
end

function add_cutting_circular_ordering_inequality!(model)
    τ = model[:τ]
    n = model[:n]
    tau = [i == j ? 0.0 : value.(model[:τ][i,j]) for i in 1:n, j in 1:n]
    w = [i == j ? 0.0 : tau[i,j] - value(model[:x][i,j,2]) for i in 1:n, j in 1:n]
    TAXA = Set(collect(2:n))
    circuit = [1]
    c_length = 0.
    while length(TAXA) > 2
        current_node = last(circuit)
        next_node = argmin(i -> w[i,current_node], TAXA)
        c_length += w[next_node, current_node]
        pop!(TAXA, next_node)
        push!(circuit, next_node)
    end
    next_node = argmin(i -> w[i,last(circuit)] + w[i,first(circuit)], TAXA)
    c_length += w[next_node, last(circuit)]
    pop!(TAXA, next_node)
    push!(circuit, next_node)
    c_length += w[first(circuit), last(circuit)]

    circuit = two_opt(w, circuit)
    c_length = circuit_length(tau, circuit)
    if c_length < 4*n-8
        # println(c_length, "<", 4n-8)
        arcs = union(Set([(circuit[i], circuit[i+1]) for i in 1:n-2]), Set([(last(circuit), first(circuit))]))
        @constraint model sum(τ[a[1], a[2]] for a in arcs) >= 4*n-8 + sum(model[:_x][i,j,2] for i in 1:n-1, j in i+1:n if (i,j) ∉ arcs)
        return true
    else
        return false
    end
end

function circuit_length(D, circuit)
    sum(D[circuit[i-1], circuit[i]] for i in 2:(length(circuit))) + D[first(circuit), last(circuit)]
end

function two_opt(D, circuit)
    new_circuit = deepcopy(circuit)
    best_length = circuit_length(D, circuit)
    while true
        l = best_length 
        for (i,j) in combinations(1:length(new_circuit), 2)
            next_j = j == length(circuit) ? 1 : j + 1
            tmp_circuit = [new_circuit[begin:i]; reverse(new_circuit[i+1:j]); new_circuit[next_j:end]]
            new_length = circuit_length(D, tmp_circuit)
            if new_length < best_length
                best_length = new_length
                new_circuit = tmp_circuit
            end          
        end
        best_length < l || break
    end
    return new_circuit
end

function solve_BME_model!(model; max_coi_cuts = 0)
    cuts = 0
    optimize!(model) 
    while cuts < max_coi_cuts
        new_cut = add_cutting_circular_ordering_inequality!(model)
        if new_cut
            cuts += 1
            optimize!(model)
        else
            break
        end
    end
end