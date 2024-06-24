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

function phylip_converter(path::String)
    rl = readlines(path)
    open(path, "w") do file
        counter = 0
        n = try 
            n = parse(Int, first(rl))
            @assert n == length(rl) - 1
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
                continue
            else
                splits = split(ln, ' ')
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
        end
        for line in rl[2:end]
            println(file, line)
        end
    end
end
