__precompile__()

# Optimal Enriched-Set Cover
module OptEnrichedSetCover

using DataFrames, Distributions, MathProgBase, JuMP

export SetMosaic, CoverParams, CoverProblem,
    CoverEnumerationParams, CoverCollection,
    nelements, ntiles, nsets, nmasks, nvars,
    nmasked, nunmasked, maskedset,
    tile, tiles, set, setsize,
    unmask, mask,
    score, # any conflicts
    optimize # conflicts with Optim.jl

using Ipopt
default_solver() = IpoptSolver(print_level=0)

include("sparse_mask_matrix.jl")
include("set_score.jl")
include("mosaic.jl")
include("masked_mosaic.jl")
include("cover_problem.jl")
include("cover_enumerator.jl")
include("parallel.jl")

end # module
