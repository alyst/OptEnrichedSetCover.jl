"""
Parameters for the `AbstractCoverProblem` (Optimal Enriched-Set Cover).
"""
struct CoverParams
    sel_prob::Float64           # prior probability to select the set, penalizes non-zero weights
    min_weight::Float64         # minimal non-zero set probability
    set_relevance_shape::Float64# how much set relevance affects set score, 0 = no effect
    setXset_factor::Float64     # how much set-set overlaps are penalized (setXset_score scale), 0 = no penalty
    setXset_shape::Float64      # how much set-set overlaps are penalized (setXset_score shape), 0 = no penalty
    maskXmask_factor::Float64   # how much activating overlapping sets in different masks is penalized
    maskXmask_shape::Float64    # how much activating overlapping sets in different masks is penalized

    function CoverParams(; sel_prob::Number=0.5, min_weight::Number = 1E-2,
                         set_relevance_shape::Number=1.0,
                         setXset_factor::Number=1.0, setXset_shape::Number=1.0,
                         maskXmask_factor::Number=1.0, maskXmask_shape::Number=1.0)
        (0.0 < sel_prob <= 1.0) || throw(ArgumentError("`set_prob` must be within (0,1] range"))
        (0.0 < min_weight <= 1.0) || throw(ArgumentError("`min_weight` must be within (0,1] range"))
        (0.0 <= set_relevance_shape) || throw(ArgumentError("`set_relevance_shape` must be ≥0"))
        (0.0 <= setXset_factor) || throw(ArgumentError("`setXset_factor` must be ≥0"))
        (0.0 <= setXset_shape) || throw(ArgumentError("`setXset_shape` must be ≥0"))
        (0.0 <= maskXmask_factor) || throw(ArgumentError("`maskXmask_factor` must be ≥0"))
        (0.0 <= maskXmask_shape) || throw(ArgumentError("`maskXmask_shape` must be ≥0"))
        new(sel_prob, min_weight, set_relevance_shape,
            setXset_factor, setXset_shape,
            maskXmask_factor, maskXmask_shape)
    end
end

"""
Linear component of an individual masked set score for the `CoverProblem`.
Doesn't take into account the overlap with the other selected sets.
"""
function msetscore_detached(masked::Number, set::Number, total_masked::Number, total::Number, relevance::Number, params::CoverParams)
    # FIXME is it just tail=:both for one set of parameters
    #= P-value for masked-vs-set overlap enriched =# res = logpvalue(masked, set, total_masked, total)*(relevance^params.set_relevance_shape) #-
    #= P-value for unmasked-vs-set overlap enriched =# #logpvalue(set - masked, set, total - total_masked, total)
    @assert !isnan(res) "masked=$masked set=$set total_masked=$total_masked total=$total relevance=$relevance res=NaN"
    return res
end

msetscore_detached(set::MaskedSet, mosaic::MaskedSetMosaic, params::CoverParams) =
    msetscore_detached(set.nmasked, set.nmasked + set.nunmasked,
                       nmasked(mosaic, set.mask), nelements(mosaic),
                       mosaic.original.set_relevances[set.set], params)

function varXvar_score(setXset::Real, iset::MaskedSet, jset::MaskedSet,
                       params::CoverParams, scale::Bool = false)
    if iset.set == jset.set
        # no penalty for the overlap with itself, even in different masks
        # (in which case it's zero to encourage the reuse of the same sets to cover different masks)
        return zero(typeof(setXset))
    elseif !isfinite(setXset)
        return -Inf # fix later with min_score
    else
        #if sXs < min_sXs
        #    min_sXs = sXs
        #end
        if iset.mask == jset.mask
            v = -((1.0-setXset)^params.setXset_shape-1.0)
            scale && (v *= params.setXset_factor)
        else
            # scale the original setXset score by maskXmask_factor if
            # the sets are from different masks,
            v = -((1.0-setXset)^(params.maskXmask_shape * params.setXset_shape)-1.0)
            scale && (v *= params.maskXmask_factor * params.setXset_factor)
        end
        return v
    end
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

"""
Result of `optimize(AbstractCoverProblem)`.
"""
struct CoverProblemResult{T, E}
    var2mset::Vector{Int}      # indices of masked sets in MaskedSetMosaic
    weights::Vector{Float64}        # weights of masked sets
    var_scores::Vector{Float64}
    total_score::T
    agg_total_score::Float64
    extra::E

    function CoverProblemResult(
            var2mset::AbstractVector{Int},
            weights::AbstractVector{Float64},
            var_scores::AbstractVector{Float64},
            total_score::T, agg_total_score::Float64, extra::E = nothing) where {T, E}
        length(var2mset) == length(weights) == length(var_scores) ||
            throw(ArgumentError("Lengths of cover result components do not match"))
        new{T, E}(var2mset, weights, var_scores,
                  total_score, agg_total_score, extra)
    end
end

function mset2var(result::CoverProblemResult, msetix::Int)
    varix = searchsortedlast(result.var2mset, msetix)
    if varix > 0 && result.var2mset[varix] == msetix
        return varix
    else
        return 0
    end
end

function selectvars(problem::AbstractCoverProblem,
                    mosaic::MaskedSetMosaic,
                    weights::AbstractVector{Float64};
                    selothermasks::Bool = true # propagate the selection to the selected sets in the other masks
)
    @assert length(weights) == nvars(problem)
    varixs = find(w -> w > problem.params.min_weight, weights) # > to make it work with min_weight=0
    if selothermasks
        maskedset2var = Dict(ms => i for (i, ms) in enumerate(problem.var2mset))
        setixs = Set(mosaic.maskedsets[problem.var2mset[i]].set for i in varixs)
        ext_varixs = Set{Int}(varixs)
        for i in setixs
            msets = mosaic.orig2masked[i]
            for ms in msets
                v = get(maskedset2var, ms, 0)
                (v > 0) && push!(ext_varixs, v)
            end
        end
        return collect(ext_varixs)
    else
        return varixs
    end
end

varXvar_mul(problem::AbstractCoverProblem, w::AbstractVector) =
    varXvar_mul!(similar(w), problem, w)

abstract type AbstractOptimizerParams{P <: AbstractCoverProblem} end;

problemtype(params::AbstractOptimizerParams{P}) where P = P
