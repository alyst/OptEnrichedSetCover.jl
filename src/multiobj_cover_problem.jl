# defines how to fold detailed score into N components
abstract type DetailedScoreFolding{N} end

# soft Fold2D score
const FoldedScore = NTuple{2, Float64}

BlackBoxOptim.fitness_scheme_type(scorefold::DetailedScoreFolding) =
    fitness_scheme_type(typeof(scorefold))

BlackBoxOptim.fitness_scheme_type(::Type{SF}) where SF <: DetailedScoreFolding{N} where N =
    ParetoFitnessScheme{N, Float64, true, MultiobjCoverProblemScoreAggregator{SF}}

struct MultiobjCoverProblemScoreAggregator{SF <: DetailedScoreFolding}
    score_folding::SF
    setXset_factor::Float64
    uncovered_factor::Float64
    covered_factor::Float64

    MultiobjCoverProblemScoreAggregator(params::CoverParams,
                                        score_folding::SF) where SF =
        new{SF}(score_folding,
                params.setXset_factor,
                params.uncovered_factor,
                params.covered_factor)
end

#=
# no folding of fitness components
struct MultiobjProblemNoFolding <: DetailedScoreFolding{4}
end
(::MultiobjProblemNoFolding)(fitness::DetailedScore) = fitness

(agg::MultiobjCoverProblemScoreAggregator{MultiobjProblemNoFolding})(score::DetailedScore) =
    score[1] + agg.setXset_factor * score[2] +
    agg.uncovered_factor * score[3] + agg.covered_factor * score[4]
=#

"""
Transforms the 4-component cover quality score into 2-component score that makes
the highly redundant solutions dominated by any less redundant ones.

See ["Cover score convolution"](@ref cover_score_convolution) section for the discussion.

# Arguments
  * `setXset_factor`: ``w_r``, same as in [`CoverParams`](@ref)
  * `uncovered_factor`: ``w_u``, same as in [`CoverParams`](@ref)
  * `covered_factor`: ``w_c``, same as in [`CoverParams`](@ref)
  * `ratio_threshold`: ``k_{\\max}``, defaults to 1
  * `shape`: ``\\alpha_k``, defaults to 0.5.
"""
struct MultiobjProblemSoftFold2d <: DetailedScoreFolding{2}
    setXset_factor::Float64
    uncovered_factor::Float64
    covered_factor::Float64
    ratio_threshold::Float64    # at what set_score/setXset_score ratio the folding kicks in
    shape::Float64              # controls folding smoothness (smaller are smoother)

    MultiobjProblemSoftFold2d(params::CoverParams,
                              ratio_threshold::Float64 = 1.0,
                              shape::Float64=0.5) =
        new(params.setXset_factor, params.uncovered_factor, params.covered_factor,
            ratio_threshold, shape)
end

function (f::MultiobjProblemSoftFold2d)(fitness::DetailedScore)
    a = fitness[1] +
        f.uncovered_factor * fitness[3] +
        f.covered_factor * fitness[4]
    b = fitness[2]
    if f.ratio_threshold !== nothing
        k = b/max(1, -a)
        (k <= f.ratio_threshold) && return (a, b)

        s = f.shape * (k - f.ratio_threshold)
        # 1/(1 + exp(1/s - s)) allows smooth transition from s = 0 (f(0)=0) to f(∞)=1
        t = inv(k * (f.setXset_factor + exp(inv(s) - s)))
        return (a + t * f.setXset_factor * b, (1.0 - t) * b)
    else
        return (a, b)
    end
end

(agg::MultiobjCoverProblemScoreAggregator{MultiobjProblemSoftFold2d})(score::FoldedScore) =
    score[1] + agg.setXset_factor * score[2]

#=
# Sums uncovered (3rd) and covered (4th) components of the score into a single score
# that to every score component mixes in the other 2 components with the score k
# The transform keeps the aggregated score intact
struct MultiobjProblemFold3d <: DetailedScoreFolding{3}
    k::Float64
    k_sXs::Float64
    k_u::Float64
    setXset_factor::Float64
    uncovered_factor::Float64
    covered_factor::Float64

    MultiobjProblemFold3d(params::CoverParams, k::Float64=0.1) =
        new(params.setXset_factor > 0.0 ? k : 0.0,
            params.setXset_factor > 0.0 ? k/params.setXset_factor : 0.0,
            params.uncovered_factor > 0.0 ? k/params.uncovered_factor : 0.0,
            params.setXset_factor,
            params.uncovered_factor,
            params.uncovered_factor > 0.0 ? params.covered_factor/params.uncovered_factor : 1.0)
end

function (f::MultiobjProblemFold3d)(fitness::DetailedScore)
    a = fitness[1]
    b = fitness[2]
    c = fitness[3] + f.covered_factor * fitness[4]
    return (0.5 * f.k * (f.setXset_factor * b + f.uncovered_factor * c) + (1-f.k) * a,
            0.5 * f.k_sXs * (a + f.uncovered_factor * c) + (1-f.k) * b,
            0.5 * f.k_u * (a + f.setXset_factor * b) + (1-f.k) * c)
end

(agg::MultiobjCoverProblemScoreAggregator{MultiobjProblemFold3d})(score::NTuple{3, Float64}) =
    score[1] + agg.setXset_factor * score[2] + agg.uncovered_factor * score[3]
=#

"""
Abstract multi-objective optimal Enriched-Set Cover problem.

See ["Method Description"](@ref method) for more details.
"""
abstract type MultiobjCoverProblem <: AbstractCoverProblem{DetailedScore} end

function score_scales(problem::MultiobjCoverProblem)
    w_min, w_max = extrema(problem.var_scores)
    sXs_min, sXs_max = length(problem.varXvar_scores) > 1 ? (Inf, -Inf) : (0.0, 0.0)
    return (max(w_max - w_min, 0.1),
            max(sXs_max - sXs_min, 0.1))
end

aggscore(w::AbstractVector{Float64}, problem::MultiobjCoverProblem) =
    aggscore(score(w, problem), problem.params)

# exclude given vars from the problem and
# adjust the remaining scores
function mask_vars(problem::MultiobjCoverProblem,
                   vars::AbstractVector{Int};
                   penalize_overlaps::Bool = true)
    varmask = fill(true, nvars(problem))
    varmask[vars] .= false
    v_scores = problem.var_scores[varmask]
    if penalize_overlaps
        # penalize overlapping sets
        pweights = fill(0.0, nvars(problem))
        pweights[vars] .= 1.0
        v_penalties = similar(pweights)
        minplus_bilinear!(take2, v_penalties,
                          problem.varXvar_scores,
                          pweights, pweights)
        v_scores .-= view(v_penalties, varmask) .* problem.params.setXset_factor
    end
    return varmask, v_scores
end

"""
Multi-objective optimal Enriched-Set Cover problem that uses MaskedSetMosaic
to define the coverage scores based on the set weights and masked elements coverage.

See ["Method Description"](@ref method) for more details.
"""
struct MaskedSetCoverProblem <: MultiobjCoverProblem
    params::CoverParams

    var2set::Vector{Int}

    var_scores::Vector{Float64}
    varXvar_scores::Matrix{Float64}

    tileXvar::SparseMaskMatrix      # tile X set(var)
    maskxtileXvar::SparseMatrixCSC{Float64} # tile-in-mask X set(var)
    nmasked_tile::Vector{Int}       # sum of masked elements for each tile in all masks (els counted for each mask separately)
    nunmasked_tile::Vector{Int}     # number of unmasked elements for each tile in the union of all masks (each el counted once)

    arrpool::ArrayPool{Float64}   # pool of vectors to reuse for evaluation FIXME move to evaluator

    function MaskedSetCoverProblem(params::CoverParams,
                        var2set::AbstractVector{Int},
                        var_scores::AbstractVector{Float64},
                        varXvar_scores::AbstractMatrix{Float64},
                        tileXvar::SparseMaskMatrix,
                        maskxtileXvar::SparseMatrixCSC,
                        nmasked_tile::AbstractVector{Int},
                        nunmasked_tile::AbstractVector{Int}
    )
        length(var2set) == length(var_scores) ==
        size(varXvar_scores, 1) == size(varXvar_scores, 2) == size(maskxtileXvar, 2) ==
        size(tileXvar, 2) ||
            throw(ArgumentError("var2set, var_scores and varXvar_scores counts do not match"))
        size(tileXvar, 1) == length(nmasked_tile) == length(nunmasked_tile) ||
            throw(ArgumentError("tileXvar and nmasked_tile rows count does not match"))
        new(params, var2set,
            var_scores, varXvar_scores,
            tileXvar, maskxtileXvar, nmasked_tile, nunmasked_tile,
            ArrayPool{Float64}())
    end
end

function MultiobjCoverProblem(mosaic::MaskedSetMosaic, params::CoverParams = CoverParams())
    v2set = var2set(mosaic)
    MaskedSetCoverProblem(params, mosaic.loc2glob_setix[v2set],
                          var_scores(mosaic, params, v2set),
                          varXvar_scores(mosaic, params, v2set, false),
                          varXtiles(mosaic, v2set, params)...)
end

ntiles(problem::MaskedSetCoverProblem) = length(problem.nmasked_tile)

function exclude_vars(problem::MaskedSetCoverProblem,
                      vars::AbstractVector{Int};
                      kwargs...)
    varmask, v_scores = mask_vars(problem, vars; kwargs...)
    return MaskedSetCoverProblem(problem.params,
                problem.var2set[varmask],
                v_scores, problem.varXvar_scores[varmask, varmask],
                problem.tileXvar[:, varmask], # FIXME can also remove unused tiles
                problem.maskxtileXvar[:, varmask], # FIXME can also remove unused tiles
                problem.nmasked_tile, problem.nunmasked_tile)
end

# the tuple of:
# - total number of uncovered masked elements (in all masks overlapping with the cover)
# - total number of covered unmasked elements (in all masks overlapping with the cover)
function miscover_score(w::AbstractVector{Float64}, problem::MaskedSetCoverProblem)
    __check_vars(w, problem)
    isempty(problem.nmasked_tile) && return (0.0, 0.0)
    # calculate the tile coverage weights
    wtile = fill!(acquire!(problem.arrpool, ntiles(problem)), 0.0)
    @inbounds for (varix, wvar) in enumerate(w)
        if wvar == 1.0
            wtile[view(problem.tileXvar, :, varix)] .= 1.0
        elseif wvar > 0.0
            for tileix in view(problem.tileXvar, :, varix)
                wtile[tileix] = max(wtile[tileix], wvar)
            end
        end
    end
    release!(problem.arrpool, wtile)
    wmaskxtile = acquire!(problem.arrpool, size(problem.maskxtileXvar, 1))
    res = sum(problem.nmasked_tile) - dot(problem.nmasked_tile, wtile),
          sum(maxplus_linear!(wmaskxtile, problem.maskxtileXvar, w))
    release!(problem.arrpool, wmaskxtile)
    return res
end

function var_scores(w::AbstractVector{Float64}, problem::MultiobjCoverProblem)
    __check_vars(w, problem)
    a = dot(problem.var_scores, w)
    if problem.params.setXset_factor == 0.0
        b = 0.0
    else
        # skip set interactions as it would be aggregated to zero
        setXset = fill!(acquire!(problem.arrpool, length(w)), 0.0)
        minplus_bilinear!(min, setXset, problem.varXvar_scores, w, w)
        b = -sum(setXset)
        release!(problem.arrpool, setXset)
    end
    return a, b
end

"""
    score(w::AbstractVector{Float64}, problem) -> NTuple{4, Float64}

Unfolded multiobjective score (fitness) of the OESC coverage.

* `w`: probabilities of the sets being covered

See ["Cover quality"](@ref cover_quality).
"""
function score(w::AbstractVector{Float64}, problem::MaskedSetCoverProblem)
    a, b = var_scores(w, problem)
    if problem.params.covered_factor == 0.0 &&
       problem.params.uncovered_factor == 0.0
        # skip miscover score calculation if it's excluded from aggregation
        c = 0.0
        d = 0.0
    else
        c, d = miscover_score(w, problem)
    end
    return (a, b, c, d)
end
