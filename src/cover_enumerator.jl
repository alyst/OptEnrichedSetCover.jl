"""
Enumerates the enriched-set covers of
the given set collection.
"""
immutable CoverEnumerator{T,S}
    mosaic::MaskedSetMosaic{T,S}    # the original masked set mosaic
    params::CoverParams             # parameters of the `CoverProblem`
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
function Base.convert(::Type{DataFrame}, covers::CoverCollection, mosaic::SetMosaic)
    masked_mosaic = mask(mosaic, covers.elmask) # FIXME a workaround, because masked_mosaic is very expensive to store in covers
    nsets = Int[sum(variant.weights.>0.0) for variant in covers.variants]
    set_ixs = vcat([find(variant.weights.>0.0)
                    for variant in covers.variants]...)
    orig_set_ixs = covers.setixs[set_ixs]
    cover_ixs = vcat([fill(cover_ix, nset) for (cover_ix, nset) in enumerate(nsets)]...)
    DataFrame(cover_ix = cover_ixs,
              set_ix = orig_set_ixs,
              set_id = mosaic.ix2set[orig_set_ixs],
              delta_score = vcat([fill(covers.variants[cover_ix].score - covers.variants[1].score, nset)
                                 for (cover_ix, nset) in enumerate(nsets)]...),
              nmasked = nmasked_perset(masked_mosaic)[set_ixs],
              nunmasked = nunmasked_perset(masked_mosaic)[set_ixs],
              weight = vcat([variant.weights[variant.weights.>0]
                              for variant in covers.variants]...),
              score = [setscore(covers, set_ixs[i], cover_ixs[i]) for i in eachindex(set_ixs)])
end

"""
Greedy enumeration of enriched-set covers.
* At each iteration a optimal enriched-set cover problem is being solved.
* The sets selected at current iteration are removed from further consideration.
* The process continues with the reduced collection until the result is an empty
  collection.

Returns `CoverCollection`.
"""
function Base.collect{T,S}(etor::CoverEnumerator{T,S}; setXset_penalty::Float64=-100.0, max_covers::Int = 0, max_set_score::Real = -10.0, max_cover_score_delta::Real = 1.0)
    cover_problem = CoverProblem(etor.mosaic, etor.params)
    res = CoverCollection(cover_problem, etor.mosaic)
    while max_covers == 0 || length(res) < max_covers
        cur_cover = optimize(cover_problem; ini_weights=rand(nsets(cover_problem)),
                             solver=IpoptSolver(print_level=0))
        used_setixs = find(cur_cover.weights .> 0)
        isempty(used_setixs) && break # no sets selected
        cover_pos = searchsortedlast(res.variants, cur_cover, by=cover->cover.score)+1
        if cover_pos > 1
            if abs(cur_cover.score - res.variants[cover_pos-1].score) <= 1E-3 &&
               maxabs(cur_cover.weights - res.variants[cover_pos-1].weights) <= 1E-3
                break # same solution, stop enumerating
            end
        end
        delta_score = 0.0
        if cover_pos > 1
            # not the best cover
            delta_score = cur_cover.score - res.variants[1].score
            if max_cover_score_delta > 0.0 && delta_score > max_cover_score_delta
                break # bad cover score, stop enumerating
            end
        end
        if isfinite(max_set_score) && (minimum(cover_problem.set_scores[used_setixs]) + delta_score > max_set_score)
            break # no set with the good score, stop enumerating
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
        scores_updated || break # the cover does not improve any score
        # save the current cover
        insert!(res.variants, cover_pos, cur_cover)
        # update pointers to the best covers for the sets
        for setix in eachindex(res.set_variantix)
            if res.set_variantix[setix] == -1
                res.set_variantix[setix] = cover_pos
            elseif res.set_variantix[setix] >= cover_pos
                # the variant has moved down
                res.set_variantix[setix] += 1
            end
        end
        (max_covers > 0 && length(res) >= max_covers) || break # collected enough covers
        # penalize selecting the same cover by penalizing every pair of sets from the cover
        for set1_ix in used_setixs, set2_ix in used_setixs
            if set1_ix != set2_ix
                cover_problem.setXset_scores[set1_ix, set2_ix] = setXset_penalty
            end
        end
        all(x::Int -> x > 0, res.set_variantix) && break # no unassigned sets left
    end
    return res
end
