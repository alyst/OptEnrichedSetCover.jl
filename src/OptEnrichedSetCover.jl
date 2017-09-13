__precompile__()

# Optimal Enriched-Set Cover
module OptEnrichedSetCover

using DataFrames, Distributions, MathProgBase, JuMP, Ipopt
using Compat

export SetMosaic, CoverParams, CoverProblem,
    CoverEnumerationParams, CoverCollection,
    nelements, ntiles, nsets, nmasks,
    nmasked, nunmasked, nmasked_pertile, nunmasked_pertile, nmasked_perset, nunmasked_perset,
    tile, tiles, set, setsize,
    unmask, mask,
    score, # any conflicts
    optimize # conflicts with Optim.jl

for script_file in ["set_score.jl",
                    "sparse_mask_matrix.jl",
                    "mosaic.jl",
                    "cover_problem.jl",
                    "cover_enumerator.jl",
                    "parallel.jl"]
    include(script_file)
end

end # module
