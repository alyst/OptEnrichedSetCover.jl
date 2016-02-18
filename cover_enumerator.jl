"""
  Enumerates the set collection coverings.
"""
immutable CoverEnumerator{T,S}
    mosaic::MaskedSetMosaic{T,S}    # the original masked set mosaic
    params::CoverParams             # parameters of the `CoverProblem`
end

"""
  The collection of masked set coverings.
"""
immutable CoverCollection
    elmask::BitVector                 # FIXME the elements mask, a workaround to avoid copying the whole mosaic upon serialization
    set_scores::Vector{Float64}
    variants::Vector{CoverProblemResult}

    function CoverCollection(mosaic::MaskedSetMosaic)
        new(mosaic.elmask,
            fill(NaN, nsets(mosaic.original)),
            Vector{CoverProblemResult}())
    end
end

function _setscore(cover::CoverProblemResult, mosaic::MaskedSetMosaic, problem::CoverProblem, setix::Int)
    if cover.weights[setix] > 0.0
        orig_setix = mosaic.setixs[setix]
        singletonsetscore(mosaic.original.set_sizes[orig_setix], nmasked_perset(mosaic)[setix],
                          nelements(mosaic), nmasked(mosaic),
                          problem.params) - log(cover.weights[setix])
    else
        Inf
    end
end

Base.length(covers::CoverCollection) = length(covers.variants)
Base.isempty(covers::CoverCollection) = isempty(covers.variants)

function Base.convert(::Type{DataFrame}, covers::CoverCollection, mosaic::SetMosaic)
    masked_mosaic = mask(mosaic, covers.elmask) # FIXME a workaround, because masked_mosaic is very expensive to store in covers
    nsets = Int[sum(variant.weights.>0.0) for variant in covers.variants]
    set_ixs = vcat([find(variant.weights.>0.0)
                    for variant in covers.variants]...)
    orig_set_ixs = masked_mosaic.setixs[set_ixs]
    DataFrame(cover_ix = vcat([fill(cover_ix, nset) for (cover_ix, nset) in enumerate(nsets)]...),
              set_ix = orig_set_ixs,
              set_id = masked_mosaic.original.ix2set[orig_set_ixs],
              delta_score = vcat([fill(covers.variants[cover_ix].score - covers.variants[1].score, nset)
                                 for (cover_ix, nset) in enumerate(nsets)]...),
              nmasked = nmasked_perset(masked_mosaic)[set_ixs],
              nunmasked = nunmasked_perset(masked_mosaic)[set_ixs],
              weight = vcat([variant.weights[variant.weights.>0]
                              for variant in covers.variants]...),
              score = covers.set_scores[orig_set_ixs])
end

"""
  Greedy enumeration of coverings.
  Sets selected at each iteration are removed from the collection.
"""
function Base.collect{T,S}(etor::CoverEnumerator{T,S}; setXset_penalty::Float64=-100.0, max_covers::Int = 0, max_set_score::Real = -10.0, max_cover_score_delta::Real = 1.0)
    cover_problem = CoverProblem(etor.mosaic, etor.params)
    res = CoverCollection(etor.mosaic)
    while max_covers == 0 || length(res) < max_covers
        cur_cover = optimize(cover_problem; ini_weights=rand(nsets(cover_problem)),
                             solver=IpoptSolver(print_level=0))
        used_setixs = find(cur_cover.weights .> 0)
        if isempty(used_setixs)
            break # no sets selected
        end
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
        if isfinite(max_set_score) && (minimum(cover_problem.set_scores[used_setixs]) > max_set_score)
            break # no set with good score, stop enumerating
        end
        scores_updated = false
        # update the set scores
        @inbounds for setix in used_setixs
            orig_setix = etor.mosaic.setixs[setix]
            # adjust the set scores by delta
            score = cover_problem.set_scores[setix] + delta_score
            if (cur_cover.weights[setix]>0.0) && isfinite(score) &&
               (isnan(res.set_scores[orig_setix]) || res.set_scores[orig_setix] > score)
               scores_updated = true
                res.set_scores[orig_setix] = score
            end
        end
        if !scores_updated
            break # the cover does not improve any score
        end
        # save current cover
        insert!(res.variants, cover_pos, cur_cover)
        if max_covers > 0 && length(res) >= max_covers
            break # collected enough covers
        end
        # penalize selective the same cover by penalizing every pair of sets from the cover
        for set1_ix in used_setixs
            for set2_ix in used_setixs
                if set1_ix != set2_ix
                    cover_problem.setXset_scores[set1_ix, set2_ix] = setXset_penalty;
                end
            end
        end
        if all(isfinite, res.set_scores)
            break # no unassigned sets left
        end
    end
    return res
end
