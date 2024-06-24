import Random.randperm
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

for n in 2 .^ [5:11;]
    D = generate_symmetric_doubly_stochastic(n)
    dirpath=joinpath("data", "RDSM$(n)")
    if !isdir(dirpath)
        mkdir(dirpath)
    end
    datapath = joinpath(dirpath,"RDSM$(n).txt")
    open(datapath, "w") do f
        println(f, string(size(D)[1]))
        for (i,r) in enumerate(eachrow(D))
            print(f,i," ")
            for e in r
                print(f, string(round(e, digits = 10) , " "))
            end
            println(f,"")    
        end
    end
    BalancedMinimumEvolution.phylip_converter(datapath)
end
