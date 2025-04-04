module BalancedMinimumEvolution

using Gurobi
const GRB_ENV_REF = Ref{Gurobi.Env}()
function __init__()
        global GRB_ENV_REF
        GRB_ENV_REF[] = Gurobi.Env()
        @eval const OPTIMIZER = () -> Gurobi.Optimizer(GRB_ENV_REF[])
        return
end

# using HiGHS
# function __init__()
#     @eval const OPTIMIZER = HiGHS.Optimizer
# end

include("graph.jl")
include("linear_program.jl")
include("neighbor_joining.jl")
include("fastme.jl")
include("beam_search.jl")
include("phylip_convert.jl")
# include("LNS/heuristic.jl")
# include("LNS/neighborhood_samplers.jl")
end # module BalancedMinimumEvolution
