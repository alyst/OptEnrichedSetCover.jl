"""
Parameters for `collect(mosaic::MaskedSetMosaic)`.
"""
struct CoverEnumerationParams
    max_covers::Int
    max_set_score::Float64
    max_cover_score_delta::Float64

    CoverEnumerationParams(;
        max_covers::Int = 0,
        max_set_score::Real = -10.0,
        max_cover_score_delta::Real = 1.0) =
        new(max_covers, max_set_score, max_cover_score_delta)
end

"""
The collection of masked set covers.
"""
struct CoverCollection
    cover_params::CoverParams
    enum_params::CoverEnumerationParams
    total_masked::Vector{Int}         # total masked elements in the mosaic
    elmasks::BitMatrix                # FIXME the elements mask, a workaround to avoid copying the whole mosaic upon serialization
    maskedsets::Vector{MaskedSet}     # sets of the MaskedSetMosaic
    base_setscores::Vector{Float64}   # base set scores
    var2cover::Vector{Int}            # best-scoring cover for the given masked set
    results::Vector{CoverProblemResult}

    function CoverCollection(problem::CoverProblem, mosaic::MaskedSetMosaic,
                             params::CoverEnumerationParams)
        # check the problem is compatible with the mosaic
        nvars(problem) == nsets(mosaic) || throw(ArgumentError("CoverProblem is not compatible to the MaskedSetMosaic: number of sets differ"))
        #nmasks(problem) == nmasks(mosaic) || throw(ArgumentError("CoverProblem is not compatible to the MaskedSetMosaic: number of masks differ"))
        new(problem.params, params, mosaic.total_masked, mosaic.elmasks, mosaic.maskedsets,
            copy(problem.set_scores),
            zeros(Int, nvars(problem)),
            Vector{CoverProblemResult}())
    end
end

function setscore(covers::CoverCollection, varix::Int, cover::CoverProblemResult, delta_score::Float64)
    if cover.weights[varix] > 0.0
        # problem set score + delta score for the best variant, where it was covered - log(set weight)
        return covers.base_setscores[varix] + delta_score - log(cover.weights[varix])
    else
        # the set is not selected
        return Inf
    end
end

setscore(covers::CoverCollection, varix::Int, coverix::Int) =
    setscore(covers, varix, covers.results[coverix],
             covers.results[coverix].total_score - covers.results[1].total_score)

function setscore(covers::CoverCollection, varix::Int)
    if covers.var2cover[varix] > 0
        # problem set score + delta score for the best variant, where it was covered - log(set weight)
        return setscore(covers, varix, covers.var2cover[varix])
    else
        # the set not selectd
        return Inf
    end
end

nmasks(covers::CoverCollection) = size(covers.elmasks, 2)
Base.length(covers::CoverCollection) = length(covers.results)
Base.isempty(covers::CoverCollection) = isempty(covers.results)

"""
Convert `covers`, a collection of the covers of `mosaic`, into a `DataFrame`.
"""
function DataFrames.DataFrame(covers::CoverCollection, mosaic::SetMosaic)
    nselsets = sum(count(x -> x > 0.0, variant.weights) for variant in covers.results)
    set_ixs = sizehint!(Vector{Int}(), nselsets)
    mask_ixs = sizehint!(Vector{Int}(), nselsets)
    cover_ixs = sizehint!(Vector{Int}(), nselsets)
    delta_scores = sizehint!(Vector{Float64}(), nselsets)
    weights = sizehint!(Vector{Float64}(), nselsets)
    scores = sizehint!(Vector{Float64}(), nselsets)
    nmasked_v = sizehint!(Vector{Int}(), nselsets)
    nunmasked_v = sizehint!(Vector{Int}(), nselsets)
    for (cover_ix, cover) in enumerate(covers.results)
        delta_score = cover.total_score - covers.results[1].total_score
        @inbounds for var_ix in eachindex(cover.weights)
            weight = cover.weights[var_ix]
            maskedset = covers.maskedsets[var_ix]
            (weight > 0.0) || continue
            push!(set_ixs, maskedset.set)
            push!(weights, weight)
            push!(mask_ixs, maskedset.mask)
            push!(cover_ixs, cover_ix)
            push!(delta_scores, delta_score)
            push!(scores, setscore(covers, var_ix, cover_ix))
            push!(nmasked_v, maskedset.nmasked)
            push!(nunmasked_v, maskedset.nunmasked)
        end
    end
    DataFrame(cover_ix = cover_ixs,
              set_ix = set_ixs,
              set_id = mosaic.ix2set[set_ixs],
              mask_ix = mask_ixs,
              delta_score = delta_scores,
              nmasked = nmasked_v,
              nunmasked = nunmasked_v,
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
    cover_coll = CoverCollection(cover_problem, mosaic, params)
    # thresholds for identifying duplicate covers
    const score_threshold = 1E-3
    const weight_threshold = 1E-3
    while true
        verbose && info("Trying to find cover #$(length(cover_coll)+1)...")
        cur_cover = optimize(cover_problem; ini_weights=rand(nvars(cover_problem)),
                             solver=default_solver())
        verbose && info("New cover found (score=$(cur_cover.total_score)), processing...")
        used_varixs = find(w -> w > 0.0, cur_cover.weights)
        if isempty(used_varixs)
            verbose && info("Cover is empty")
            break
        end
        cover_pos = searchsortedlast(cover_coll.results, cur_cover, by=cover->cover.total_score)+1
        if cover_pos > 1
            cover = cover_coll.results[cover_pos-1]
            if abs(cur_cover.total_score - cover.total_score) <= score_threshold &&
               all(i -> (@inbounds return abs(cur_cover.weights[i] - cover.weights[i]) <= weight_threshold),
                   eachindex(cur_cover.weights))
                verbose && info("Duplicate solution")
                break
            end
        end
        delta_score = 0.0
        if cover_pos > 1
            # not the best cover
            delta_score = cur_cover.total_score - cover_coll.results[1].total_score
            if params.max_cover_score_delta > 0.0 && delta_score > params.max_cover_score_delta
                verbose && info("Cover score_delta=$(delta_score) above threshold")
                break
            end
        end
        if isfinite(params.max_set_score) &&
           all(i -> (@inbounds return all(x -> x + delta_score > params.max_set_score, view(cover_problem.set_scores, i, :))),
               used_varixs)
            verbose && info("All set scores below $(max_set_score)")
            break
        end
        scores_updated = false
        # update the best set scores
        @inbounds for varix in used_varixs
            # adjust the set scores by delta
            new_score = setscore(cover_coll, varix, cur_cover, delta_score)
            cur_score = setscore(cover_coll, varix)
            if isfinite(new_score) && (!isfinite(cur_score) || cur_score > new_score)
                scores_updated = true
                cover_coll.var2cover[varix] = -1 # mark for setting to the cover_pos
            end
        end
        if !scores_updated
            verbose && info("No global set scores improvement")
            break
        end
        # save the current cover
        insert!(cover_coll.results, cover_pos, cur_cover)
        verbose && info("Cover saved")
        # update pointers to the best covers for the sets
        for sXmix in eachindex(cover_coll.var2cover)
            if cover_coll.var2cover[sXmix] == -1
                cover_coll.var2cover[sXmix] = cover_pos
            elseif cover_coll.var2cover[sXmix] >= cover_pos
                # the variant has moved down
                cover_coll.var2cover[sXmix] += 1
            end
        end
        if params.max_covers > 0 && length(cover_coll) >= params.max_covers
            verbose && info("Maximal number of covers collected")
            break
        end
        # penalize selecting the sets from the current cover again
        cover_problem.set_scores[used_varixs] = -log(cover_params.sel_prob)
        if all(x::Int -> x > 0, cover_coll.var2cover)
            verbose && info("All sets assigned to covers")
            break
        end
    end
    verbose && info("$(length(cover_coll)) cover(s) collected")
    return cover_coll
end
