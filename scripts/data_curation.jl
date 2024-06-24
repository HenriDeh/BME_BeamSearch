using BalancedMinimumEvolution

for dataset in [["01-Primates12", "02-M17", "03-M18", "04-SeedPlants500", "05-M43", "06-M62", "07-RbcL55", "08-Rana64", "09-M82"];["100_rdpii_F81", "100_rdpii_F84", "100_rdpii_K2P", "100_rdpii_JC69","200_rdpii_F81", "200_rdpii_F84", "200_rdpii_K2P", "200_rdpii_JC69","300_zilla_F81", "300_zilla_F84", "300_zilla_K2P", "300_zilla_JC69"]; ["RDSM32"]; ["RDSM64"];["RDSM128"];["RDSM256"];["RDSM512"];["RDSM1024"];["RDSM2048"]]
    path = joinpath("data",dataset,dataset*".txt")
    BalancedMinimumEvolution.to_phylip(path)
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

run(`find . -name "*.nwk*" `)
