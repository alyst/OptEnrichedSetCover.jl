__precompile__()

# Optimal Enriched-Set Cover
module OptEnrichedSetCover

using DataFrames, Distributions, MathProgBase, JuMP, BlackBoxOptim

export SetMosaic, CoverParams,
    AbstractCoverProblem, QuadraticCoverProblem, MultiobjectiveCoverProblem,
    AbstractOptimizerParams, QuadraticOptimizerParams, MultiobjectiveOptimizerParams,
    CoverEnumerationParams, CoverCollection, CoverProblemResult,
    nelements, ntiles, nsets, nmaskedsets, nmasks, nvars,
    nmasked, nunmasked, maskedset,
    tile, tiles, set, setsize,
    unmask, mask,
    aggscore, score, # any conflicts
    optimize # conflicts with Optim.jl

using Ipopt
default_quadratic_solver() = IpoptSolver(print_level=0)
#using Mosek
#default_quadratic_solver() = MosekSolver()
#using Gurobi
#default_quadratic_solver() = GurobiSolver(OutputFlag=0)
#using CPLEX
#default_quadratic_solver() = CplexSolver()

import BlackBoxOptim: OptimizationProblem

include("sparse_mask_matrix.jl")
include("set_score.jl")
include("mosaic.jl")
include("masked_mosaic.jl")
include("cover_problem.jl")
include("quadratic_cover_problem.jl")
include("multiobj_cover_problem.jl")
include("cover_enumerator.jl")
include("parallel.jl")

end # module
