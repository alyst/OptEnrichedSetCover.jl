"""
Parameters for `collect(mosaic::MaskedSetMosaic)`.
"""
struct CoverEnumerationParams
    max_covers::Int
    max_set_score::Union{Float64, Void}
    max_cover_score_delta::Union{Float64, Void}

    CoverEnumerationParams(;
        max_covers::Int = 0,
        max_set_score::Union{Real, Void} = -10.0,
        max_cover_score_delta::Union{Real, Void} = nothing) =
        new(max_covers, max_set_score, max_cover_score_delta)
end

"""
The collection of masked set covers.
"""
struct CoverCollection{T}
    cover_params::CoverParams
    enum_params::CoverEnumerationParams
    total_masked::Vector{Int}         # total masked elements in the mosaic
    elmasks::BitMatrix                # FIXME the elements mask, a workaround to avoid copying the whole mosaic upon serialization
    maskedsets::Vector{MaskedSet}     # sets of the MaskedSetMosaic
    base_msetscores::Vector{Float64}  # base set scores
    var2cover::Vector{Int}            # best-scoring cover for the given masked set
    results::Vector{CoverProblemResult{T}}

    function CoverCollection(problem::AbstractCoverProblem{T}, mosaic::MaskedSetMosaic,
                             params::CoverEnumerationParams) where T
        # check the problem is compatible with the mosaic
        nvars(problem) == nsets(mosaic) || throw(ArgumentError("CoverProblem is not compatible to the MaskedSetMosaic: number of sets differ"))
        #nmasks(problem) == nmasks(mosaic) || throw(ArgumentError("CoverProblem is not compatible to the MaskedSetMosaic: number of masks differ"))
        new{T}(problem.params, params, mosaic.total_masked, mosaic.elmasks, mosaic.maskedsets,
               copy(problem.var_scores),
               zeros(Int, nvars(problem)),
               Vector{CoverProblemResult{T}}())
    end
end

function msetscore(covers::CoverCollection, varix::Int, cover::CoverProblemResult)
    if cover.weights[varix] > 0.0
        # problem set score + delta score for the best variant, where it was covered - log(set weight)
        return covers.base_msetscores[varix] - log(cover.weights[varix])
    else
        # the set is not selected
        return Inf
    end
end

msetscore(covers::CoverCollection, varix::Int, coverix::Int) =
    msetscore(covers, varix, covers.results[coverix])

function msetscore(covers::CoverCollection, varix::Int)
    if covers.var2cover[varix] > 0
        # problem set score + delta score for the best variant, where it was covered - log(set weight)
        return msetscore(covers, varix, covers.var2cover[varix])
    else
        # the set is not covered
        return NaN
    end
end

nmasks(covers::CoverCollection) = size(covers.elmasks, 2)
Base.length(covers::CoverCollection) = length(covers.results)
Base.isempty(covers::CoverCollection) = isempty(covers.results)

function OptimizerParams(problem_type::Symbol; kwargs...)
    if problem_type == :quadratic
        return QuadraticOptimizerParams(kwargs...)
    elseif problem_type == :multiobjective
        return MultiobjectiveOptimizerParams(kwargs...)
    else
        throw(ArgumentError("Unsupported cover problem type $problem_type"))
    end
end

"""
Greedy enumeration of enriched-set covers.
* At each iteration an optimal enriched-set cover problem is being solved.
* The sets selected at current iteration are removed from further consideration.
* The process continues with the reduced collection until the result is an empty
  collection or the last cover is much worse than the first one.

Returns `CoverCollection`.
"""
Base.collect(mosaic::MaskedSetMosaic,
             cover_params::CoverParams=CoverParams(),
             enum_params::CoverEnumerationParams=CoverEnumerationParams();
             problem_type::Symbol=:quadratic,
             verbose::Bool=false,
             optargs...) =
    collect(mosaic, cover_params, enum_params, OptimizerParams(problem_type, optargs...), verbose)

function Base.collect(mosaic::MaskedSetMosaic,
                      cover_params::CoverParams,
                      enum_params::CoverEnumerationParams,
                      opt_params::AbstractOptimizerParams,
                      verbose::Bool)
    cover_problem = problemtype(opt_params)(mosaic, cover_params)
    collect!(CoverCollection(cover_problem, mosaic, enum_params),
             cover_problem, mosaic, enum_params, opt_params, verbose)
end

function collect!(cover_coll::CoverCollection,
                  cover_problem::AbstractCoverProblem,
                  mosaic::MaskedSetMosaic,
                  params::CoverEnumerationParams,
                  opt_params::AbstractOptimizerParams,
                  verbose::Bool=false
)
    verbose && info("Starting covers enumeration...")
    # thresholds for identifying duplicate covers
    const score_threshold = 1E-3
    const weight_threshold = 1E-3
    while true
        verbose && info("Trying to find cover #$(length(cover_coll)+1)...")
        cur_cover = optimize(cover_problem, opt_params)
        verbose && info("New cover found (score=$(cur_cover.total_score), aggscore=$(cur_cover.agg_total_score)), processing...")
        used_varixs = find(w -> w > 0.0, cur_cover.weights)
        if isempty(used_varixs)
            verbose && info("Cover is empty")
            break
        end
        cover_pos = searchsortedlast(cover_coll.results, cur_cover, by=cover->cover.agg_total_score)+1
        if cover_pos > 1
            cover = cover_coll.results[cover_pos-1]
            if abs(cur_cover.agg_total_score - cover.agg_total_score) <= score_threshold &&
               all(i -> (@inbounds return abs(cur_cover.weights[i] - cover.weights[i]) <= weight_threshold),
                   eachindex(cur_cover.weights))
                verbose && info("Duplicate solution")
                break
            end
        end
        if cover_pos > 1 && isa(params.max_cover_score_delta, Float64)
            # not the best cover
            cover_score_delta = cur_cover.agg_total_score - cover_coll.results[1].agg_total_score
            if cover_score_delta > params.max_cover_score_delta
                verbose && info("Cover score_delta=$(cover_score_delta) above threshold")
                break
            end
        end
        if isa(params.max_set_score, Float64)
            max_set_score = params.max_set_score::Float64
            if all(i -> (@inbounds return all(x -> x > max_set_score, view(cover_problem.var_scores, i, :))),
                    used_varixs)
                verbose && info("All set scores below $(max_set_score)")
                break
            end
        end
        scores_updated = false
        # update the best set scores
        @inbounds for varix in used_varixs
            # adjust the set scores by delta
            new_score = msetscore(cover_coll, varix, cur_cover)
            cur_score = msetscore(cover_coll, varix)
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
        verbose && info("Cover collected")
        # update pointers to the best covers for the sets
        for i in eachindex(cover_coll.var2cover)
            if cover_coll.var2cover[i] == -1
                cover_coll.var2cover[i] = cover_pos
            elseif cover_coll.var2cover[i] >= cover_pos
                # the variant has moved down
                cover_coll.var2cover[i] += 1
            end
        end
        if params.max_covers > 0 && length(cover_coll) >= params.max_covers
            verbose && info("Maximal number of covers collected")
            break
        end
        if all(x -> x > 0, cover_coll.var2cover)
            verbose && info("All sets assigned to covers")
            break
        end
        # exclude current solution
        penalize_solution!(cover_problem, mosaic, cur_cover.weights)
    end
    verbose && info("$(length(cover_coll)) cover(s) collected")
    return cover_coll
end

"""
Convert `covers`, a collection of the covers of `mosaic`, into a `DataFrame`.
"""
function DataFrames.DataFrame(covers::CoverCollection, mosaic::SetMosaic; report::Symbol=:covered)
    if report == :covered
        return report_covered(covers, mosaic)
    elseif report == :matrix
        return report_matrix(covers, mosaic)
    else
        throw(ArgumentError("Unknown report mode '$report'"))
    end
end

# DataFrame report for the covered sets (weight > 0)
function report_covered(covers::CoverCollection, mosaic::SetMosaic)
    nselsets = sum(count(x -> x > 0.0, variant.weights) for variant in covers.results)
    set_ixs = sizehint!(Vector{Int}(), nselsets)
    mask_ixs = sizehint!(Vector{Int}(), nselsets)
    cover_ixs = sizehint!(Vector{Int}(), nselsets)
    cover_scores = sizehint!(Vector{Float64}(), nselsets)
    weights = sizehint!(Vector{Float64}(), nselsets)
    relevances = sizehint!(Vector{Float64}(), nselsets)
    cv_scores = sizehint!(Vector{Float64}(), nselsets)
    de_scores = sizehint!(Vector{Float64}(), nselsets)
    nmasked_v = sizehint!(Vector{Int}(), nselsets)
    nunmasked_v = sizehint!(Vector{Int}(), nselsets)
    for (cover_ix, cover) in enumerate(covers.results)
        @inbounds for var_ix in eachindex(cover.weights)
            weight = cover.weights[var_ix]
            maskedset = covers.maskedsets[var_ix]
            (weight > 0.0) || continue
            push!(set_ixs, maskedset.set)
            push!(weights, weight)
            push!(relevances, mosaic.set_relevances[maskedset.set])
            push!(mask_ixs, maskedset.mask)
            push!(cover_ixs, cover_ix)
            push!(cover_scores, cover.agg_total_score)
            push!(cv_scores, msetscore(covers, var_ix, cover))
            push!(de_scores, msetscore_detached(maskedset.nmasked, setsize(maskedset),
                                                covers.total_masked[maskedset.mask], nelements(mosaic),
                                                1.0 #= ignore relevance =#, covers.cover_params))
            push!(nmasked_v, maskedset.nmasked)
            push!(nunmasked_v, maskedset.nunmasked)
        end
    end
    DataFrame(cover_ix = cover_ixs,
              set_ix = set_ixs,
              set_id = mosaic.ix2set[set_ixs],
              mask_ix = mask_ixs,
              cover_score = cover_scores,
              nmasked = nmasked_v,
              nunmasked = nunmasked_v,
              relevance = relevances,
              mset_cover_weight = weights,
              mset_score_covered = cv_scores,
              mset_score_detached = de_scores)
end

# DataFrame report for each mask X each set,
# where all set that are covered in at least one mask are considered
# only the best covers used
function report_matrix(covers::CoverCollection, mosaic::SetMosaic)
    # collect all sets that are covered in at least one mask
    selsets = Set{Int}()
    for (varix, coverix) in enumerate(covers.var2cover)
        if coverix > 0
            push!(selsets, covers.maskedsets[varix].set)
        end
    end
    nmasks = size(covers.elmasks, 2)
    sets_v = sort!(collect(selsets))
    setsizes_v = setsize.(mosaic, sets_v)
    set2index = Dict(zip(sets_v, 1:length(sets_v)))
    nmasked_mtx = nmasked_perset(mosaic, covers.elmasks, set2index)
    coverix_mtx = zeros(Int, size(nmasked_mtx))
    weights_mtx = zeros(Float64, size(nmasked_mtx))
    cover_scores_mtx = zeros(Float64, size(nmasked_mtx))
    cv_scores_mtx = fill(NaN, size(nmasked_mtx))
    @inbounds for (varix, coverix) in enumerate(covers.var2cover)
        coverix > 0 || continue
        mset = covers.maskedsets[varix]
        setix = set2index[mset.set]
        i = sub2ind(size(coverix_mtx), setix, mset.mask)
        coverix_mtx[i] = coverix
        weights_mtx[i] = covers.results[coverix].weights[varix]
        cover_scores_mtx[i] = covers.results[coverix].agg_total_score
        cv_scores_mtx[i] = msetscore(covers, varix, coverix)
    end
    de_scores_mtx = zeros(Float64, size(nmasked_mtx))
    @inbounds for j in 1:nmasks
        for i in eachindex(sets_v)
            setix = sets_v[i]
            de_scores_mtx[i, j] = msetscore_detached(nmasked_mtx[i, j], setsizes_v[i],
                                                     covers.total_masked[j], nelements(mosaic),
                                                     1.0 #= ignore relevance =#, covers.cover_params)
        end
    end
    DataFrame(cover_ix = vec(coverix_mtx),
              set_ix = repeat(sets_v, outer=[nmasks]),
              set_id = repeat(mosaic.ix2set[sets_v], outer=[nmasks]),
              mask_ix = repeat(collect(1:nmasks), inner=[length(sets_v)]),
              cover_score = vec(cover_scores_mtx),
              nmasked = vec(nmasked_mtx),
              nunmasked = repeat(setsizes_v, outer=[nmasks]) .- vec(nmasked_mtx),
              relevance = repeat(mosaic.set_relevances[sets_v], outer=[nmasks]),
              cover_weight = vec(weights_mtx),
              mset_score_covered = vec(cv_scores_mtx),
              mset_score_detached = vec(de_scores_mtx))
end
