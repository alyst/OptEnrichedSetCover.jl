"""
Parameters for the `AbstractCoverProblem` (Optimal Enriched-Set Cover).
"""
struct CoverParams
    sel_prob::Float64           # prior probability to select the set, penalizes non-zero weights
    min_weight::Float64         # minimal non-zero set probability
    set_relevance_shape::Float64# how much set relevance affects set score, 0 = no effect
    setXset_factor::Float64     # how much set-set overlaps are penalized (setXset_score scale), 0 = no penalty
    setXset_shape::Float64      # how much set-set overlaps are penalized (setXset_score shape), 0 = no penalty
    uncovered_factor::Float64   # how much masked uncovered elements penalize the score
    covered_factor::Float64     # how much unmasked covered elements penalize the score

    function CoverParams(; sel_prob::Number=0.5, min_weight::Number = 1E-2,
                         set_relevance_shape::Number=1.0,
                         setXset_factor::Number=1.0, setXset_shape::Number=1.0,
                         uncovered_factor::Number=0.1, covered_factor::Number=0.025)
        (0.0 < sel_prob <= 1.0) || throw(ArgumentError("`set_prob` must be within (0,1] range"))
        (0.0 < min_weight <= 1.0) || throw(ArgumentError("`min_weight` must be within (0,1] range"))
        (0.0 <= set_relevance_shape) || throw(ArgumentError("`set_relevance_shape` must be ≥0"))
        (0.0 <= setXset_factor) || throw(ArgumentError("`setXset_factor` must be ≥0"))
        (0.0 <= setXset_shape) || throw(ArgumentError("`setXset_shape` must be ≥0"))
        (0.0 <= uncovered_factor) || throw(ArgumentError("`uncovered_factor` must be ≥0"))
        (0.0 <= covered_factor) || throw(ArgumentError("`covered_factor` must be ≥0"))
        new(sel_prob, min_weight, set_relevance_shape,
            setXset_factor, setXset_shape,
            uncovered_factor, covered_factor)
    end
end

"""
Linear component of an individual masked set score for the `CoverProblem`.
Doesn't take into account the overlap with the other selected sets.
"""
function overlap_score(masked::Number, set::Number, total_masked::Number, total::Number,
                       relevance::Number, params::CoverParams)
    # FIXME is it just tail=:both for one set of parameters
    #= P-value for masked-vs-set overlap enriched =# res = logpvalue(masked, set, total_masked, total)*(relevance^params.set_relevance_shape) #-
    #= P-value for unmasked-vs-set overlap enriched =# #logpvalue(set - masked, set, total - total_masked, total)
    @assert !isnan(res) "masked=$masked set=$set total_masked=$total_masked total=$total relevance=$relevance res=NaN"
    return res
end

overlap_score(set::MaskedSet, mosaic::MaskedSetMosaic, params::CoverParams) =
    overlap_score(set.nmasked, set.nmasked + set.nunmasked,
                  nmasked(mosaic, set.mask), nelements(mosaic),
                  mosaic.original.set_relevances[set.set], params)

var2set(mosaic::MaskedSetMosaic) = sort!(collect(keys(mosaic.orig2masked)))

# extracts
# 1) var-to-tile mapping (only relevant tiles are considered)
# 2) the total number of masked elements in all masks for each tile
# 3) the total number of unmasked elements in all active masks (overlapping with the tile) for each tile
function tilemaskXvar(mosaic::MaskedSetMosaic)
    setixs = var2set(mosaic)
    nm = nmasks(mosaic)
    used_tilemasks = Set{Tuple{Int, Int}}()
    @inbounds for (varix, setix) in enumerate(setixs)
        set_tiles = view(mosaic.original.tileXset, :, setix)
        for msetix in mosaic.orig2masked[setix]
            maskix = mosaic.maskedsets[msetix].mask
            for tileix in set_tiles
                push!(used_tilemasks, (tileix, maskix))
            end
        end
    end
    used_tilemask_v = sort!(collect(used_tilemasks)) # sort lexicographically
    nmasked_pertilemask = fill(0, length(used_tilemask_v))
    nunmasked_pertilemask = fill(0, length(used_tilemask_v))
    last_tileix = 0
    tile_size = 0
    tile_elms = view([], 1:0) # empty range for type stability
    for (ix, (tileix, maskix)) in enumerate(used_tilemask_v)
        if tileix != last_tileix
            tile_elms = view(mosaic.original.elmXtile, :, tileix)
            tile_size = length(tile_elms)
            last_tileix = tileix
        end
        nmasked = sum(view(mosaic.elmasks, tile_elms, maskix))
        nmasked_pertilemask[ix] = nmasked
        nunmasked_pertilemask[ix] = tile_size - nmasked
    end
    tileXvar = mosaic.original.tileXset[first.(used_tilemask_v), setixs]
    if isempty(tileXvar)
        return tileXvar, nmasked_pertilemask, nunmasked_pertilemask
    end

    # find tiles that refer to the identical set of variables
    varXtile = transpose(tileXvar)
    vars2tileixs = Dict{Vector{Int}, Vector{Int}}()
    for tileix in 1:size(varXtile, 2)
        varixs = varXtile[:, tileix]
        tileixs = get!(() -> Vector{Int}(), vars2tileixs, varixs)
        push!(tileixs, tileix)
    end
    if length(vars2tileixs) == size(varXtile, 2) # no duplicates
        return tileXvar, nmasked_pertilemask, nunmasked_pertilemask
    end
    tileXvar = tileXvar[[first(tileixs) for tileixs in values(vars2tileixs)], :] # take the 1st of duplicated tiles
    # aggregate n[un]masked over duplicated tiles
    nmasked_pertilegroup = [sum(i -> nmasked_pertilemask[i], tileixs)
                            for tileixs in values(vars2tileixs)]
    nunmasked_pertilegroup = [sum(i -> nunmasked_pertilemask[i], tileixs)
                            for tileixs in values(vars2tileixs)]
    return tileXvar, nmasked_pertilegroup, nunmasked_pertilegroup
end

# var score is:
# the sum of overlap scores with all masked sets
# minus the log(var selection probability)
# plus nets*log(nsets), where nsets arr all masked sets overlapping with var
function var_scores(mosaic::MaskedSetMosaic, var2set::AbstractVector{Int}, params::CoverParams)
    # calculate the sum of scores of given set in each mask
    sel_penalty = log(params.sel_prob)
    v_scores = Vector{Float64}(undef, length(var2set))
    for (varix, setix) in enumerate(var2set)
        scoresum = 0.0
        msetixs = mosaic.orig2masked[setix]
        for msetix in msetixs
            scoresum += overlap_score(mosaic.maskedsets[msetix], mosaic, params)
        end
        v_scores[varix] = scoresum - sel_penalty + length(msetixs)*log(length(msetixs))
    end
    return v_scores
end

function varXvar_score(setXset::Real, params::CoverParams, scale::Bool = false)
    vXv = 1.0-(1.0-setXset)^params.setXset_shape
    scale && (vXv *= params.setXset_factor)
    return vXv
end

function varXvar_scores(mosaic::MaskedSetMosaic, var2set::AbstractVector{Int},
                        params::CoverParams, scale::Bool = false)
    vXv_scores = varXvar_score.(mosaic.original.setXset_scores[var2set, var2set], Ref(params), scale)
    for i in eachindex(var2set) # don't consider self-intersections
        @inbounds vXv_scores[i, i] = zero(eltype(vXv_scores))
    end
    vXv_min = Inf # minimum finite varXvar_scores element
    @inbounds for vXv in vXv_scores
        if isfinite(vXv) && vXv < vXv_min
            vXv_min = vXv
        end
    end
    # replace infinite varXvar score with the minimal finite setXset score
    vXv_subst = varXvar_score(1.25*vXv_min, params, true)
    @inbounds for i in eachindex(vXv_scores)
        if !isfinite(vXv_scores[i])
            v1, v2 = ind2sub(size(vXv_scores), i)
            @warn("var[$v1]×var[$v2] score is $(vXv_scores[i])")
            if vXv_scores[i] < 0.0
                vXv_scores[i] = vXv_subst
            end
        end
    end
    return vXv_scores
end

"""
Optimal Enriched-Set Cover problem -- choose the sets from the collection `𝒞` to cover
the masked(selected) elements `M`.
The optimal sets cover `C = {c₁, c₂, ..., cₙ} ⊂ 𝒞` has to deliver 3 goals:
* be relevant (i.e. minimize the P-values of `M` and `cᵢ` sets overlap)
* be minimal (i.e. minimize the number of sets in `C`)
* be non-redundant (i.e. minimize the P-values of the pairwise non-overlap of `C` sets with each other).

Fuzzy set selection is possible -- each set is assigned a weight from `[0, 1]` range.
"""
abstract type AbstractCoverProblem{T} end

"""
Total number of problem variables.
"""
nvars(problem::AbstractCoverProblem) = length(problem.var_scores)

# clamps weights to 0..1 range and returns sum(|w - fixed_w|)
function fix_cover_weights!(weights::Matrix{Float64})
    pen = 0.0
    @inbounds for i in eachindex(weights)
        w = uncov_probs[i]
        w_new = clamp(w, 0.0, 1.0)
        pen += abs2(w_new - w)
        weights[i] = w_new
    end
    return pen
end

__check_vars(w::AbstractVector{Float64}, problem::AbstractCoverProblem) =
    length(w) == nvars(problem) ||
        throw(DimensionMismatch("Incorrect length of parameters vector: $(length(w)) ($(nvars(problem)) expected)"))

"""
Result of `optimize(AbstractCoverProblem)`.
"""
struct CoverProblemResult{T}
    var2set::Vector{Int}        # indices of sets in the original SetMosaic
    weights::Vector{Float64}    # solution: weights of the sets
    var_scores::Vector{Float64} # scores of the sets
    total_score::T
    agg_total_score::Float64
    extra

    function CoverProblemResult(
            var2set::AbstractVector{Int},
            weights::AbstractVector{Float64},
            var_scores::AbstractVector{Float64},
            total_score::T, agg_total_score::Float64, extra = nothing) where T
        length(var2set) == length(weights) == length(var_scores) ||
            throw(DimensionMismatch("Lengths of cover result components do not match"))
        new{T}(var2set, weights, var_scores,
               total_score, agg_total_score, extra)
    end
end

function set2var(result::CoverProblemResult, setix::Int)
    varix = searchsortedlast(result.var2set, setix)
    if varix > 0 && result.var2set[varix] == setix
        return varix
    else
        return 0
    end
end

function selectvars(problem::AbstractCoverProblem,
                    mosaic::MaskedSetMosaic,
                    weights::AbstractVector{Float64}
)
    @assert length(weights) == nvars(problem)
    return findall(w -> w > problem.params.min_weight, weights) # > to make it work with min_weight=0
end

abstract type AbstractOptimizerParams{P <: AbstractCoverProblem} end;

problemtype(params::AbstractOptimizerParams{P}) where P = P
