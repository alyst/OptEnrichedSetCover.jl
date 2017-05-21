__precompile__(true)

# Optimal Enriched-Set Cover
module OptEnrichedSetCover

using DataFrames, Distributions, MathProgBase, JuMP, Ipopt

export SetMosaic, CoverParams, CoverProblem,
    CoverEnumerationParams, CoverCollection,
    nelements, ntiles, nsets,
    nmasked, nunmasked, nmasked_pertile, nunmasked_pertile, nmasked_perset, nunmasked_perset,
    tile, tiles, set,
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
