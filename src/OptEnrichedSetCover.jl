__precompile__()

# Optimal Enriched-Set Cover
module OptEnrichedSetCover

using DataFrames, Distributions, MathProgBase, JuMP, Ipopt

export SetMosaic, CoverParams, CoverProblem,
    CoverEnumerationParams, CoverCollection,
    nelements, ntiles, nsets, nmasks,
    nmasked, nunmasked, nmasked_pertile, nunmasked_pertile, nmasked_perset, nunmasked_perset,
    tile, tiles, set, setsize,
    unmask, mask,
    score, # any conflicts
    optimize # conflicts with Optim.jl

include("set_score.jl")
include("sparse_mask_matrix.jl")
include("mosaic.jl")
include("masked_mosaic.jl")
include("cover_problem.jl")
include("cover_enumerator.jl")
include("parallel.jl")

end # module
