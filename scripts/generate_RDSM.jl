include("../src/phylip_convert.jl")

function generate_symmetric_doubly_stochastic(n; tol=1e-9, max_iter=1000)
    A = rand(n, n)
    A = (A + A') / 2
    for i in 1:n
        A[i, i] = 0
    end
    # Normalize the matrix to be doubly stochastic using Sinkhorn-Knopp algorithm
    for _ in 1:max_iter
        row_sums = sum(A, dims=2)
        A = A ./ row_sums 
        A = (A + A') / 2

        col_sums = sum(A, dims=1)
        A = A ./ col_sums
        A = (A + A') / 2
        # Check for convergence
        if all(abs.(row_sums .- 1) .< tol) && all(abs.(col_sums .- 1) .< tol)
            break
        end
    end

    return A
end

for n in 10:5:60
    D = generate_symmetric_doubly_stochastic(n)
    dirpath=joinpath("data", "RDSM$(n)")
    if !isdir(dirpath)
        mkdir(dirpath)
    end
    datapath = joinpath(dirpath,"RDSM$(n).txt")
    D_to_txt(datapath, D)
end
