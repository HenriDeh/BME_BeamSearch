using Printf
ndecdigits(v) = ceil(Int, -log10(min(eps(v), eps(typeof(v)))))
function expand_float_string(s) 
    if !occursin("e", s) && length(s) < 23
        return s
    else
        d = parse(Float64,s)
        s = @sprintf "%1.*f" min(ndecdigits(d),19) d
        return s
    end
end

function to_phylip(path::String)
    rl = readlines(path)
    file = open(path, "w") #do file
    try 
        counter = 0
        header = false
        n = try 
            n = parse(Int, first(rl))
            @assert n == length(rl) - 1
            header = true
            n
        catch e
            length(rl)
        end
        println(file, n)
        for (i,ln) in enumerate(rl)
            if ln == ""
                counter = 0
                continue
            end
            if counter == 0
                counter += 1
                if header
                    continue
                end
            end
            splits = split(ln, ' ')
            filter!(!=(""), splits)
            if !startswith(ln, string(counter)) || length(splits) > n
                if length(splits) > n
                    splits = splits[(1 + length(splits) - n):end]
                end
                rl[i] = "$counter "*join(map(expand_float_string, splits), " ")
            else
                rl[i] = join(map(expand_float_string, splits), " ")
            end
            counter += 1
        end
        for line in rl[end-n+1:end]
            println(file, line)
        end
    finally
        close(file)
    end
end


function D_to_txt(datapath, D::Matrix)
    open(datapath, "w") do f
        println(f, string(size(D,1)))
        for (i,r) in enumerate(eachrow(D))
            print(f,i," ")
            for e in r
                print(f, string(round(e, digits = 10) , " "))
            end
            println(f,"")    
        end
    end
    to_phylip(datapath)
end
