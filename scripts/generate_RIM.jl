include("../src/phylip_convert.jl")

function generate_symmetric_integer_matrix(N::Int)
    mat = zeros(Int, N, N)
    for i in 1:N
        for j in i+1:N
            mat[i, j] = rand(1:10)
            mat[j, i] = mat[i, j]
        end
    end
    return mat
end

for n in 10:5:60
    D = generate_symmetric_integer_matrix(n)
    dirpath = joinpath("data", "RIM$(n)")
    if !isdir(dirpath)
        mkdir(dirpath)
    end
    datapath = joinpath(dirpath, "RIM$(n).txt")
    D_to_txt(datapath, D)
end