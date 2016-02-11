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
immutable CoverCollection{T,S}
    mosaic::MaskedSetMosaic{T,S}
    setscores::Vector{Float64}
    variants::Vector{CoverProblemResult}

    function CoverCollection(mosaic::MaskedSetMosaic{T,S})
        new(mosaic,
            fill(NaN, nsets(mosaic.original)),
            Vector{CoverProblemResult}())
    end
end

function _setscore(covers::CoverCollection, cover::CoverProblemResult, problem::CoverProblem, setix::Int)
    if cover.weights[setix] > 0.0
        orig_setix = covers.mosaic.setixs[setix]
        delta_score = isempty(covers.variants) ? 0.0 : cover.score - covers.variants[1].score
        setscore(covers.mosaic.original.set_sizes[orig_setix], nmasked_perset(covers.mosaic)[setix],
                nelements(covers.mosaic), nmasked(covers.mosaic),
                problem.params) - log(cover.weights[setix]) + delta_score
    else
        Inf
    end
end

Base.length(covers::CoverCollection) = length(covers.variants)
Base.isempty(covers::CoverCollection) = isempty(covers.variants)

function Base.convert(::Type{DataFrame}, covers::CoverCollection)
    nsets = Int[sum(variant.weights.>0.0) for variant in covers.variants]
    set_ixs = vcat([find(variant.weights.>0.0)
                    for variant in covers.variants]...)
    orig_set_ixs = covers.mosaic.setixs[set_ixs]
    DataFrame(cover_ix = vcat([fill(cover_ix, nset) for (cover_ix, nset) in enumerate(nsets)]...),
              set_ix = orig_set_ixs,
              set_id = covers.mosaic.original.ix2set[orig_set_ixs],
              delta_score = vcat([fill(covers.variants[cover_ix].score - covers.variants[1].score, nset)
                                 for (cover_ix, nset) in enumerate(nsets)]...),
              nmasked = nmasked_perset(covers.mosaic)[set_ixs],
              nunmasked = nunmasked_perset(covers.mosaic)[set_ixs],
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
    res = CoverCollection{T,S}(etor.mosaic)
    best_cover_score = NaN
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
        set_scores = [_setscore(res, cur_cover, cur_problem, i) for i in eachindex(cur_cover.weights)]
        cover_score = score(orig_problem, orig_weights)
        if isempty(res)
            best_cover_score = cover_score
        elseif max_cover_score_delta > 0.0 && cover_score - best_cover_score > max_cover_score_delta
            break # bad cover score
        end
        if isfinite(max_set_score) && (minimum(set_scores) > max_set_score)
            break # no set with good score
        end
        # save current cover
        push!(res.variants, CoverProblemResult(orig_weights, cover_score))
        # update the set scores
        for (setix, score) in enumerate(set_scores)
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
