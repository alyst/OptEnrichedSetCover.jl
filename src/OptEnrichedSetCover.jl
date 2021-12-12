# Optimal Enriched-Set Cover
module OptEnrichedSetCover

using Requires, DataFrames, Distributions, BlackBoxOptim
using LinearAlgebra, SparseArrays, Distributed, StatsBase
using Printf: @printf

export SetMosaic, CoverParams,
    AbstractCoverProblem, MultiobjCoverProblem,
    AbstractOptimizerParams, MultiobjOptimizerParams,
    CoverEnumerationParams, CoverCollection, CoverProblemResult,
    nelements, ntiles, nsets, nmaskedsets, nmasks, nexperiments, nvars, nsolutions,
    nmasked, nunmasked, maskedset, overlap,
    tile, tiles, set, setsize,
    setidtype, experimentidtype, weighttype,
    originalmosaic, mask,
    set_relevance,
    logpvalue, aggscore, score, # any conflicts
    best_aggscore, best_index, best_varweights, varweights,
    optimize # conflicts with Optim.jl

import BlackBoxOptim: OptimizationProblem

include("array_pool.jl")
include("sparse_mask_matrix.jl")
include("set_score.jl")
include("mosaic.jl")
include("weighted_mosaic.jl")
include("masked_mosaic.jl")
include("cover_problem.jl")
include("multiobj_cover_problem_utils.jl")
include("multiobj_cover_problem.jl")
include("multiobj_cover_problem_bbo.jl")
include("cover_enumerator.jl")
include("parallel.jl")

end # module
