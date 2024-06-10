# Generating random matrices
for n in [20,50,100,200,400,800,1600]
    M = rand(n,n)
    dirpath = joinpath("data","Rand$n")
    if !isdir(dirpath)
        mkdir(dirpath)
    end
    path = joinpath(dirpath, "Rand$n.txt")
    open(path, "w") do f
        println(f, n)
        for (i,l) in enumerate(eachrow(M))
            print(f, i, " ")
            println(f, join(round.(l, digits = 7), " "))
        end
    end
    phylip_converter(path)
end

for dataset in ["01-Primates12","02-M17","03-M18","04-SeedPlants500","05-M43","06-M62","07-RbcL55","08-Rana64","09-M82"]
    path = joinpath("data",dataset,dataset*".txt")
    phylip_converter(path)
end

#Clean datafiles WARNING: this erases all newick files. Only datasets are kept.
for (root, dirs, files) in walkdir("data")
    for exp_dir in dirs
        for (root2, exps, datafiles) in walkdir(joinpath(root,exp_dir))
            for exp in exps
                rm(joinpath(root2,exp), recursive = true)
            end
            for file in datafiles
                if file ∉ ["01-Primates12","02-M17","03-M18","04-SeedPlants500","05-M43","06-M62","07-RbcL55","08-Rana64","09-M82"] .* ".txt"
                    rm(joinpath(root2,file))
                end
            end
        end
    end
end