import Random.randperm
include("../src/phylip_convert.jl")
function randomDoublyStochasticMatrix(n, num_perm = n^2)
    M = zeros(n,n)
    α = rand(num_perm)
    α = α / sum(α)

    for i=1:num_perm
        perm = randperm(n)
        for j=1:n
            M[perm[j],j] += α[i]
        end
    end

    return M
end

for n in 2 .^ [5:11;]
    for i in 1:50
        D = randomDoublyStochasticMatrix(n, Int(ceil(n*log(n))))
        D = (D .+ D')./2
        dirpath=joinpath("data", "RDSM$(n)_$i")
        if !isdir(dirpath)
            mkdir(dirpath)
        end
        datapath = joinpath(dirpath,"RDSM$(n)_$i.txt")
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
        phylip_converter(datapath)
    end
end