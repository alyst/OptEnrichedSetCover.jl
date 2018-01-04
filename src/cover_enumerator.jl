"""
Parameters for `collect(mosaic::MaskedSetMosaic)`.
"""
struct CoverEnumerationParams
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
struct CoverCollection
    total_masked::Vector{Int}         # total masked elements in the mosaic
    elmasks::BitMatrix                # FIXME the elements mask, a workaround to avoid copying the whole mosaic upon serialization
    setixs::Vector{Int}               # vector of the set indices in the original mosaic
    base_setscores::Matrix{Float64}   # base set scores
    var2cover::Matrix{Int}          # best-scoring cover for the given set
    results::Vector{CoverProblemResult}

    CoverCollection(empty::Void = nothing) =
        new(Vector{Int}(), BitVector(), Vector{Int}(), Vector{Float64}(),
            Vector{Int}(), Vector{CoverProblemResult}())

    function CoverCollection(problem::CoverProblem, mosaic::MaskedSetMosaic)
        nsets(problem) == nsets(mosaic) || throw(ArgumentError("CoverProblem is not compatible to the MaskedSetMosaic: number of sets differ"))
        nmasks(problem) == nmasks(mosaic) || throw(ArgumentError("CoverProblem is not compatible to the MaskedSetMosaic: number of masks differ"))
        # FIXME check the problem is compatible with the mosaic
        new(copy(mosaic.total_masked), mosaic.elmasks, mosaic.setixs,
            problem.set_scores,
            zeros(Int, nsets(problem), nmasks(problem)),
            Vector{CoverProblemResult}())
    end
end

function setscore(covers::CoverCollection, setix::Int, maskix::Int, cover::CoverProblemResult, delta_score::Float64)
    if cover.weights[setix, maskix] > 0.0
        # problem set score + delta score for the best cover, where it was covered - log(set weight)
        return covers.base_setscores[setix, maskix] + delta_score - log(cover.weights[setix, maskix])
    else
        # the set is not selectd
        return Inf
    end
end

setscore(covers::CoverCollection, setix::Int, maskix::Int, coverix::Int) =
    setscore(covers, setix, maskix, covers.results[coverix],
             covers.results[coverix].score - covers.results[1].score)

function setscore(covers::CoverCollection, setix::Int, maskix::Int)
    if covers.var2cover[setix, maskix] > 0
        # problem set score + delta score for the best cover, where it was covered - log(set weight)
        return setscore(covers, setix, maskix, covers.var2cover[setix])
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
    # restore masked mosaic
    # FIXME a workaround, because masked_mosaic is very expensive to store in covers
    masked_mosaic = MaskedSetMosaic(mosaic, covers.elmasks, covers.setixs)
    nsets = Int[count(x -> x > 0.0, view(cover.weights, :, maskix)) for maskix in 1:nmasks(covers), cover in covers.results]
    nsetsum = sum(nsets)
    set_ixs = sizehint!(Vector{Int}(), nsetsum)
    mask_ixs = sizehint!(Vector{Int}(), nsetsum)
    cover_ixs = sizehint!(Vector{Int}(), nsetsum)
    delta_scores = sizehint!(Vector{Float64}(), nsetsum)
    weights = sizehint!(Vector{Float64}(), nsetsum)
    scores = sizehint!(Vector{Float64}(), nsetsum)
    nmasked_v = sizehint!(Vector{Int}(), nsetsum)
    nunmasked_v = sizehint!(Vector{Int}(), nsetsum)
    for (cover_ix, cover) in enumerate(covers.results)
        delta_score = cover.score - covers.results[1].score
        @inbounds for set_ix in 1:size(cover.weights, 1), mask_ix in 1:size(cover.weights, 2)
            weight = cover.weights[set_ix, mask_ix]
            (weight > 0.0) || continue
            push!(set_ixs, set_ix)
            push!(weights, weight)
            push!(mask_ixs, mask_ix)
            push!(cover_ixs, cover_ix)
            push!(delta_scores, delta_score)
            push!(scores, setscore(covers, set_ix, mask_ix, cover_ix))
            push!(nmasked_v, masked_mosaic.nmasked_perset[set_ix, mask_ix])
            push!(nunmasked_v, masked_mosaic.nunmasked_perset[set_ix, mask_ix])
        end
    end
    orig_set_ixs = covers.setixs[set_ixs]
    DataFrame(cover_ix = cover_ixs,
              set_ix = orig_set_ixs,
              set_id = mosaic.ix2set[orig_set_ixs],
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
    res = CoverCollection(cover_problem, mosaic)
    # thresholds for identifying duplicate covers
    const score_threshold = 1E-3
    const weight_threshold = 1E-3
    while true
        verbose && info("Trying to find cover #$(length(res)+1)...")
        cur_cover = optimize(cover_problem; ini_weights=rand(nsets(cover_problem)),
                             solver=IpoptSolver(print_level=0))
        verbose && info("New cover found (score=$(cur_cover.score)), processing...")
        used_setixs = filter(setix -> any(w -> w > 0.0, view(cur_cover.weights, setix, :))::Bool,
                             1:size(cur_cover.weights, 1))
        if isempty(used_setixs)
            verbose && info("Cover is empty")
            break
        end
        cover_pos = searchsortedlast(res.results, cur_cover, by=cover->cover.score)+1
        if cover_pos > 1
            cover = res.results[cover_pos-1]
            if abs(cur_cover.score - cover.score) <= score_threshold &&
               all(i -> (@inbounds return abs(cur_cover.weights[i] - cover.weights[i]) <= weight_threshold),
                   eachindex(cur_cover.weights))
                verbose && info("Duplicate solution")
                break
            end
        end
        delta_score = 0.0
        if cover_pos > 1
            # not the best cover
            delta_score = cur_cover.score - res.results[1].score
            if params.max_cover_score_delta > 0.0 && delta_score > params.max_cover_score_delta
                verbose && info("Cover score_delta=$(delta_score) above threshold")
                break
            end
        end
        if isfinite(params.max_set_score) &&
           all(i -> (@inbounds return all(x -> x + delta_score > params.max_set_score, view(cover_problem.set_scores, i, :))),
               used_setixs)
            verbose && info("All set scores below $(max_set_score)")
            break
        end
        scores_updated = false
        # update the set scores
        @inbounds for setix in used_setixs, maskix in 1:nmasks(cover_problem)
            # adjust the set scores by delta
            new_score = setscore(res, setix, maskix, cur_cover, delta_score)
            cur_score = setscore(res, setix, maskix)
            if isfinite(new_score) && (!isfinite(cur_score) || cur_score > new_score)
                scores_updated = true
                res.var2cover[setix, maskix] = -1 # mark for setting to the cover_pos
            end
        end
        if !scores_updated
            verbose && info("No global set scores improvement")
            break
        end
        # save the current cover
        insert!(res.results, cover_pos, cur_cover)
        verbose && info("Cover saved")
        # update pointers to the best covers for the sets
        for sXmix in eachindex(res.var2cover)
            if res.var2cover[sXmix] == -1
                res.var2cover[sXmix] = cover_pos
            elseif res.var2cover[sXmix] >= cover_pos
                # the cover has moved down
                res.var2cover[sXmix] += 1
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
        if all(x::Int -> x > 0, res.var2cover)
            verbose && info("All sets assigned to covers")
            break
        end
    end
    verbose && info("$(length(res)) cover(s) collected")
    return res
end
