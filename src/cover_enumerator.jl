"""
Parameters for `collect(mosaic::MaskedSetMosaic)`.
"""
immutable CoverEnumerationParams
    max_covers::Int
    max_set_score::Float64
    max_cover_score_delta::Float64
    setXset_penalty::Float64

    CoverEnumerationParams(;
        max_covers::Int = 0,
        max_set_score::Real = -10.0,
        max_cover_score_delta::Real = 1.0,
        setXset_penalty::Float64=-100.0) =
        new(max_covers, max_set_score, max_cover_score_delta, setXset_penalty)
end

"""
The collection of masked set covers.
"""
immutable CoverCollection
    elmask::BitVector                 # FIXME the elements mask, a workaround to avoid copying the whole mosaic upon serialization
    setixs::Vector{Int}               # vector of the set indices in the original mosaic
    base_setscores::Vector{Float64}   # base set scores
    set_variantix::Vector{Int}        # best-scoring cover for the given set
    variants::Vector{CoverProblemResult}

    CoverCollection(empty::Void = nothing) =
        new(BitVector(), Vector{Int}(), Vector{Float64}(),
            Vector{Int}(), Vector{CoverProblemResult}())

    function CoverCollection(problem::CoverProblem, mosaic::MaskedSetMosaic)
        nsets(problem) == nsets(mosaic) || throw(ArgumentError("CoverProblem is not compatible to the MaskedSetMosaic: number of sets differ"))
        # FIXME check the problem is compatible with the mosaic
        new(mosaic.elmask, mosaic.setixs,
            problem.set_scores,
            zeros(Int, nsets(problem)),
            Vector{CoverProblemResult}())
    end
end

function setscore(covers::CoverCollection, setix::Int, cover::CoverProblemResult, delta_score::Float64)
    if cover.weights[setix] > 0.0
        # problem set score + delta score for the best variant, where it was covered - log(set weight)
        return covers.base_setscores[setix] + delta_score - log(cover.weights[setix])
    else
        # the set is not selectd
        return Inf
    end
end

setscore(covers::CoverCollection, setix::Int, variantix::Int) = setscore(covers, setix, covers.variants[variantix],
                                                                         covers.variants[variantix].score - covers.variants[1].score)

function setscore(covers::CoverCollection, setix::Int)
    if covers.set_variantix[setix] > 0
        # problem set score + delta score for the best variant, where it was covered - log(set weight)
        return setscore(covers, setix, covers.set_variantix[setix])
    else
        # the set not selectd
        return Inf
    end
end

Base.length(covers::CoverCollection) = length(covers.variants)
Base.isempty(covers::CoverCollection) = isempty(covers.variants)

"""
Convert `covers`, a collection of the covers of `mosaic`, into a `DataFrame`.
"""
function DataFrames.DataFrame(covers::CoverCollection, mosaic::SetMosaic)
    # restore masked mosaic
    # FIXME a workaround, because masked_mosaic is very expensive to store in covers
    masked_mosaic = MaskedSetMosaic(mosaic, covers.elmask, covers.setixs)
    nsets = Int[count(x -> x > 0.0, variant.weights) for variant in covers.variants]
    nsetsum = sum(nsets)
    set_ixs = sizehint!(Vector{Int}(), nsetsum)
    cover_ixs = sizehint!(Vector{Int}(), nsetsum)
    delta_scores = sizehint!(Vector{Float64}(), nsetsum)
    weights = sizehint!(Vector{Float64}(), nsetsum)
    scores = sizehint!(Vector{Float64}(), nsetsum)
    for (variant_ix, variant) in enumerate(covers.variants)
        delta_score = variant.score - covers.variants[1].score
        for (set_ix, weight) in enumerate(variant.weights)
            (weight > 0.0) || continue
            push!(set_ixs, set_ix)
            push!(weights, weight)
            push!(cover_ixs, variant_ix)
            push!(delta_scores, delta_score)
            push!(scores, setscore(covers, set_ix, variant_ix))
        end
    end
    orig_set_ixs = covers.setixs[set_ixs]
    DataFrame(cover_ix = cover_ixs,
              set_ix = orig_set_ixs,
              set_id = mosaic.ix2set[orig_set_ixs],
              delta_score = delta_scores,
              nmasked = nmasked_perset(masked_mosaic)[set_ixs],
              nunmasked = nunmasked_perset(masked_mosaic)[set_ixs],
              weight = weights,
              score = scores)
end

"""
Greedy enumeration of enriched-set covers.
* At each iteration an optimal enriched-set cover problem is being solved.
* The sets selected at current iteration are removed from further consideration.
* The process continues with the reduced collection until the result is an empty
  collection or the last cover is much worse than the first one.

Returns `CoverCollection`.
"""
function Base.collect(mosaic::MaskedSetMosaic,
                      cover_params::CoverParams=CoverParams(),
                      params::CoverEnumerationParams=CoverEnumerationParams();
                      verbose::Bool=false
)
    verbose && info("Starting covers enumeration...")
    cover_problem = CoverProblem(mosaic, cover_params)
    res = CoverCollection(cover_problem, mosaic)
    while true
        verbose && info("Trying to find cover #$(length(res)+1)...")
        cur_cover = optimize(cover_problem; ini_weights=rand(nsets(cover_problem)),
                             solver=IpoptSolver(print_level=0))
        verbose && info("New cover found (score=$(cur_cover.score)), processing...")
        used_setixs = find(w -> w > 0.0, cur_cover.weights)
        if isempty(used_setixs)
            verbose && info("Cover is empty")
            break
        end
        cover_pos = searchsortedlast(res.variants, cur_cover, by=cover->cover.score)+1
        if cover_pos > 1
            variant = res.variants[cover_pos-1]
            if abs(cur_cover.score - variant.score) <= 1E-3 &&
               maxabs(cur_cover.weights - variant.weights) <= 1E-3
                verbose && info("Duplicate solution")
                break
            end
        end
        delta_score = 0.0
        if cover_pos > 1
            # not the best cover
            delta_score = cur_cover.score - res.variants[1].score
            if params.max_cover_score_delta > 0.0 && delta_score > params.max_cover_score_delta
                verbose && info("Cover score_delta=$(delta_score) above threshold")
                break
            end
        end
        if isfinite(params.max_set_score) && (minimum(cover_problem.set_scores[used_setixs]) + delta_score > params.max_set_score)
            verbose && info("All set scores below $(max_set_score)")
            break
        end
        scores_updated = false
        # update the set scores
        @inbounds for setix in used_setixs
            # adjust the set scores by delta
            new_score = setscore(res, setix, cur_cover, delta_score)
            cur_score = setscore(res, setix)
            if isfinite(new_score) && (!isfinite(cur_score) || cur_score > new_score)
                scores_updated = true
                res.set_variantix[setix] = -1 # mark for setting to the cover_pos
            end
        end
        if !scores_updated
            verbose && info("No global set scores improvement")
            break
        end
        # save the current cover
        insert!(res.variants, cover_pos, cur_cover)
        verbose && info("Cover saved")
        # update pointers to the best covers for the sets
        for setix in eachindex(res.set_variantix)
            if res.set_variantix[setix] == -1
                res.set_variantix[setix] = cover_pos
            elseif res.set_variantix[setix] >= cover_pos
                # the variant has moved down
                res.set_variantix[setix] += 1
            end
        end
        if params.max_covers > 0 && length(res) >= params.max_covers
            verbose && info("Maximal number of covers collected")
            break
        end
        # penalize selecting the same cover by penalizing every pair of sets from the cover
        for set1_ix in used_setixs, set2_ix in used_setixs
            if set1_ix != set2_ix
                cover_problem.setXset_scores[set1_ix, set2_ix] = params.setXset_penalty
            end
        end
        if all(x::Int -> x > 0, res.set_variantix)
            verbose && info("All sets assigned to covers")
            break
        end
    end
    verbose && info("$(length(res)) cover(s) collected")
    return res
end
