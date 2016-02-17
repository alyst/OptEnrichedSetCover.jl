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
    setscores::Vector{Float64}
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
              score = covers.setscores[orig_set_ixs])
end

"""
  Greedy enumeration of coverings.
  Sets selected at each iteration are removed from the collection.
"""
function Base.collect{T,S}(etor::CoverEnumerator{T,S}; max_covers::Int = 0, max_set_score::Real = -10.0, max_cover_score_delta::Real = 1.0)
    orig_problem = CoverProblem(etor.mosaic, etor.params)
    cur_mosaic = copy(etor.mosaic)
    res = CoverCollection(etor.mosaic)
    while max_covers == 0 || length(res) < max_covers
        cur_problem = CoverProblem(cur_mosaic, etor.params)
        cur_cover = optimize(cur_problem; ini_weights=rand(nsets(cur_problem)),
                             solver=IpoptSolver(print_level=0))
        if sum(cur_cover.weights) == 0.0
            break # no sets selected
        end
        # get the weights of sets in the original etor.mosaic problem
        orig_weights = fill(0.0, nsets(etor.mosaic))
        orig_weights[Int[findfirst(etor.mosaic.setixs, setix) for setix in cur_mosaic.setixs]] = cur_cover.weights
        set_scores = [_setscore(cur_cover, cur_mosaic, cur_problem, i) for i in eachindex(cur_cover.weights)]
        # transform cur_cover for the same set collection as in the original problem
        orig_cover = CoverProblemResult(orig_weights, score(orig_problem, orig_weights))
        cover_pos = searchsortedlast(res.variants, orig_cover, by=cover->cover.score)+1
        if cover_pos > 1
            # not the best cover
            delta_score = orig_cover.score - res.variants[1].score
            if max_cover_score_delta > 0.0 && delta_score > max_cover_score_delta
                break # bad cover score, stop enumerating
            end
            # adjust the set scores by delta
            @inbounds for i in 1:length(set_scores)
                set_scores[i] += delta_score
            end
        end
        if isfinite(max_set_score) && (minimum(set_scores) > max_set_score)
            break # no set with good score, stop enumerating
        end
        # save current cover
        insert!(res.variants, cover_pos, orig_cover)
        # update the set scores
        @inbounds for (setix, score) in enumerate(set_scores)
            if isfinite(score)
                res.setscores[cur_mosaic.setixs[setix]] = score
            end
        end
        if max_covers > 0 && length(res) >= max_covers
            break # collected enough covers
        end
        # remove selected sets from the current collection
        filter!(cur_mosaic, cur_cover.weights .== 0)
        if nsets(cur_mosaic) == 0
            break # no sets left
        end
    end
    return res
end
