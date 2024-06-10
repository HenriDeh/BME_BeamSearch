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
        for (i,ln) in enumerate(rl)
            if ln == ""
                counter = 0
                continue
            end
            if counter == 0
                counter += 1
                continue
            else
                if !startswith(ln, string(counter))
                    rl[i] = "$counter "*join(map(expand_float_string, split(ln, ' ')), " ")
                else
                    rl[i] = join(map(expand_float_string, split(ln, ' ')), " ")
                end
                counter += 1
            end
        end
        for line in rl
            println(file, line)
        end
    end
end
