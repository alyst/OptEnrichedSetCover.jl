"""
Parameters for the `AbstractCoverProblem` (Optimal Enriched-Set Cover).
"""
struct CoverParams
    sel_prob::Float64           # prior probability to select the set, penalizes non-zero weights
    min_weight::Float64         # minimal non-zero set probability
    set_relevance_shape::Float64# how much set relevance affects set score, 0 = no effect
    set_relevance_min::Float64  # if shaped relevance is below, it's set to set_relevance_min
    setXset_factor::Float64     # how much set-set overlaps are penalized (setXset_score scale), 0 = no penalty
    uncovered_factor::Float64   # how much masked uncovered elements penalize the score
    covered_factor::Float64     # how much unmasked covered elements penalize the score

    function CoverParams(; sel_prob::Number=0.5, min_weight::Number = 1E-2,
                         set_relevance_shape::Number=0.5,
                         set_relevance_min::Number=0.5,
                         setXset_factor::Number=1.0,
                         uncovered_factor::Number=0.1, covered_factor::Number=0.025)
        (0.0 < sel_prob <= 1.0) || throw(ArgumentError("`set_prob` must be within (0,1] range"))
        (0.0 < min_weight <= 1.0) || throw(ArgumentError("`min_weight` must be within (0,1] range"))
        (0.0 <= set_relevance_shape) || throw(ArgumentError("`set_relevance_shape` must be ≥0"))
        (0.0 <= set_relevance_min <= 1) || throw(ArgumentError("`set_relevance_min` must be within [0, 1] range"))
        (0.0 <= setXset_factor) || throw(ArgumentError("`setXset_factor` must be ≥0"))
        (0.0 <= uncovered_factor) || throw(ArgumentError("`uncovered_factor` must be ≥0"))
        (0.0 <= covered_factor) || throw(ArgumentError("`covered_factor` must be ≥0"))
        new(sel_prob, min_weight, set_relevance_shape, set_relevance_min,
            setXset_factor, uncovered_factor, covered_factor)
    end
end

"""
Linear component of an individual masked set score for the `CoverProblem`.
Doesn't take into account the overlap with the other selected sets.
"""
function overlap_score(masked::Number, set::Number, total_masked::Number, total::Number,
                       relevance::Number, params::CoverParams)
    # FIXME is it just tail=:both for one set of parameters
    #= P-value for masked-vs-set overlap enriched =# res = logpvalue(masked, set, total_masked, total)*max(relevance^params.set_relevance_shape, params.set_relevance_min) #-
    #= P-value for unmasked-vs-set overlap enriched =# #logpvalue(set - masked, set, total - total_masked, total)
    @assert !isnan(res) "masked=$masked set=$set total_masked=$total_masked total=$total relevance=$relevance res=NaN"
    return res
end

overlap_score(molap::MaskOverlap, setix::Int, mosaic::MaskedSetMosaic, params::CoverParams) =
    overlap_score(molap.nmasked, molap.nmasked + molap.nunmasked,
                  nmasked(mosaic, molap.mask), nelements(mosaic),
                  mosaic.original.set_relevances[setix], params)

var2set(mosaic::MaskedSetMosaic) = sort!(collect(keys(mosaic.set2masks)))

# extracts
# 1) var-to-tile mapping (only relevant tiles are considered)
# 2) the total number of masked elements in all masks for each tile
# 3) the total number of unmasked elements in all active masks (overlapping with the tile) for each tile
function tilemaskXvar(mosaic::MaskedSetMosaic, var2setix::AbstractVector{Int} = var2set(mosaic))
    # find relevant tiles and
    # sum masked elements per tile in all masks and unmasked elements per tile in union mask
    # note: in masks there could be element not covered by any set/tile, so these counts refer to coverable elements
    tileixs = Vector{Int}()
    tile2nmasked = Vector{Int}()
    tile2nunmasked = Vector{Int}()
    @inbounds for (varix, setix) in enumerate(var2setix)
        set_tiles = view(mosaic.original.tileXset, :, setix)
        for tileix in set_tiles
            tilepos = searchsortedfirst(tileixs, tileix)
            if tilepos > length(tileixs) || tileixs[tilepos] != tileix # do this once per tile
                tile_elms = view(mosaic.original.elmXtile, :, tileix)
                insert!(tileixs, tilepos, tileix)
                insert!(tile2nunmasked, tilepos, length(tile_elms) - sum(view(mosaic.elunionmask, tile_elms)))
                nmasked = 0
                for molap in mosaic.set2masks[setix]
                    nmasked += sum(view(mosaic.elmasks, tile_elms, molap.mask))
                end
                insert!(tile2nmasked, tilepos, nmasked)
            end
        end
    end

    tileXvar = mosaic.original.tileXset[tileixs, var2setix]
    if !isempty(tileXvar)
        # find tiles that refer to the identical set of vars and merge their stats
        varXtile = transpose(tileXvar)
        vars2tiles = Dict{Vector{Int}, Vector{Int}}()
        for tileix in axes(varXtile, 2)
            varixs = varXtile[:, tileix]
            tileixs = get!(() -> Vector{Int}(), vars2tiles, varixs)
            push!(tileixs, tileix)
        end
        if length(vars2tiles) < size(tileXvar, 1) # there are tiles to merge
            tileXvar = tileXvar[[first(tileixs) for tileixs in values(vars2tiles)], :] # take the 1st of merged tiles
            # aggregate n[un]masked over duplicated tiles
            tile2nmasked = [sum(i -> tile2nmasked[i], tileixs) for tileixs in values(vars2tiles)]
            tile2nunmasked = [sum(i -> tile2nunmasked[i], tileixs) for tileixs in values(vars2tiles)]
        end
    end
    return tileXvar, tile2nmasked, tile2nunmasked
end

# var score is:
# the sum of overlap scores with all masked sets
# minus the log(var selection probability)
# plus nets*log(nsets), where nsets arr all masked sets overlapping with var
function var_scores(mosaic::MaskedSetMosaic, var2set::AbstractVector{Int}, params::CoverParams)
    # calculate the sum of scores of given set in each mask
    sel_penalty = -log(params.sel_prob)
    v_scores = Vector{Float64}(undef, length(var2set))
    @inbounds for (varix, setix) in enumerate(var2set)
        scoresum = 0.0
        for molap in mosaic.set2masks[setix]
            scoresum = overlap_score(molap, setix, mosaic, params)
        end
        v_scores[varix] = scoresum + sel_penalty
    end
    return v_scores
end

varXvar_score(setXset::Real, params::CoverParams, scale::Bool = false) =
    ifelse(scale, setXset * params.setXset_factor, setXset)

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
    vXv_cartixs = CartesianIndices(vXv_scores)
    @inbounds for i in eachindex(vXv_scores)
        if !isfinite(vXv_scores[i])
            v1, v2 = Tuple(vXv_cartixs[i])
            @warn "var[$v1]×var[$v2] score is $(vXv_scores[i])"
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

function selectvars(problem::AbstractCoverProblem,
                    mosaic::MaskedSetMosaic,
                    weights::AbstractVector{Float64})
    @assert length(weights) == nvars(problem)
    return findall(w -> w > problem.params.min_weight, weights) # > to make it work with min_weight=0
end

abstract type AbstractOptimizerParams{P <: AbstractCoverProblem} end;

problemtype(params::AbstractOptimizerParams{P}) where P = P
