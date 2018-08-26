"""
Parameters for `collect(mosaic::MaskedSetMosaic)`.
"""
struct CoverEnumerationParams
    max_covers::Int
    max_set_score::Union{Float64, Nothing}
    max_cover_score_delta::Union{Float64, Nothing}

    CoverEnumerationParams(;
        max_covers::Int = 0,
        max_set_score::Union{Real, Nothing} = -10.0,
        max_cover_score_delta::Union{Real, Nothing} = nothing) =
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
    setixs::Vector{Int}               # sets referenced by the MaskedSetMosaic
    base_selscores::Vector{Float64}   # base selected set scores
    sel2cover::Vector{Int}            # best-scoring cover for the given selected set
    results::Vector{CoverProblemResult{T}}

    function CoverCollection(problem::AbstractCoverProblem{T}, mosaic::MaskedSetMosaic,
                             params::CoverEnumerationParams) where T
        # check the problem is compatible with the mosaic
        nvars(problem) == nsets(mosaic) ||
                throw(ArgumentError("CoverProblem is not compatible to the MaskedSetMosaic: number of vars and sets differ"))
        #nmasks(problem) == nmasks(mosaic) || throw(ArgumentError("CoverProblem is not compatible to the MaskedSetMosaic: number of masks differ"))
        new{T}(problem.params, params, mosaic.total_masked, mosaic.elmasks,
               problem.var2set, copy(problem.var_scores), zeros(Int, nvars(problem)),
               Vector{CoverProblemResult{T}}())
    end
end

function set2sel(covers::CoverCollection, setix::Int, selix::Int=0)
    if selix == 0
        return searchsortedlast(covers.setixs, setix)
    elseif covers.setixs[selix] != setix
        throw(ArgumentError("Incorrect selix hint ($selix)"))
    else
        return selix
    end
end

function setscore(covers::CoverCollection, setix::Int,
                  cover::CoverProblemResult, selix::Int=0, varix::Int=0)
    selix = set2sel(covers, setix, selix)
    (selix == 0) && return NaN # the set is not a part of the cover collection
    if varix == 0
        varix = set2var(cover, setix)
        (varix == 0) && return NaN # the set is not a part of the cover problem
    elseif cover.var2set[varix] != setix
        throw(ArgumentError("Incorrect varix hint ($varix)"))
    end
    if cover.weights[varix] > 0.0
        # problem set score + delta score for the best variant, where it was covered - log(set weight)
        return covers.base_selscores[selix] - log(cover.weights[varix])
    else
        # the set is not selected
        return Inf
    end
end

function setscore(covers::CoverCollection, setix::Int, selix::Int=0)
    selix = set2sel(covers, setix, selix)
    (selix == 0) && return NaN # the set is not a part of the cover collection
    if covers.sel2cover[selix] > 0
        # problem set score + delta score for the best variant, where it was covered - log(set weight)
        return setscore(covers, setix, covers.sel2cover[selix], selix)
    else
        # the set is not covered
        return NaN
    end
end

nmasks(covers::CoverCollection) = size(covers.elmasks, 2)
Base.length(covers::CoverCollection) = length(covers.results)
Base.isempty(covers::CoverCollection) = isempty(covers.results)

function OptimizerParams(problem_type::Symbol, params::ParamsDict)
    if problem_type == :quadratic
        return QuadraticOptimizerParams(;params...)
    elseif problem_type == :multiobjective
        return MultiobjectiveOptimizerParams(;params...)
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
    collect(mosaic, cover_params, enum_params,
            OptimizerParams(problem_type, ParamsDict(optargs)),
            verbose)

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
    verbose && @info("Starting covers enumeration...")
    # thresholds for identifying duplicate covers
    score_threshold = 1E-3
    weight_threshold = 1E-3
    while true
        verbose && @info("Trying to find cover #$(length(cover_coll)+1)...")
        cur_cover = optimize(cover_problem, opt_params)
        verbose && @info("New cover found (score=$(cur_cover.total_score), aggscore=$(cur_cover.agg_total_score)), processing...")
        used_varixs = findall(w -> w > cover_problem.params.min_weight, cur_cover.weights)
        if isempty(used_varixs)
            verbose && @info("Cover is empty")
            break
        end
        # store cover in the list maintaining sorting order
        cover_pos = searchsortedlast(cover_coll.results, cur_cover, by=cover->cover.agg_total_score)+1
        if cover_pos > 1
            cover = cover_coll.results[cover_pos-1]
            if abs(cur_cover.agg_total_score - cover.agg_total_score) <= score_threshold &&
               cur_cover.var2set == cover.var2set #= should never happen =# &&
               all(i -> (@inbounds return abs(cur_cover.weights[i] - cover.weights[i]) <= weight_threshold),
                   eachindex(cur_cover.weights))
                verbose && @info("Duplicate solution")
                break
            end
        end
        if cover_pos > 1 && isa(params.max_cover_score_delta, Float64)
            # not the best cover
            cover_score_delta = cur_cover.agg_total_score - cover_coll.results[1].agg_total_score
            if cover_score_delta > params.max_cover_score_delta
                verbose && @info("Cover score_delta=$(cover_score_delta) above threshold")
                break
            end
        end
        if isa(params.max_set_score, Float64)
            max_set_score = params.max_set_score::Float64
            if all(i -> (@inbounds return all(x -> x > max_set_score, view(cover_problem.var_scores, i, :))),
                    used_varixs)
                verbose && @info("All set scores below $(max_set_score)")
                break
            end
        end
        scores_updated = false
        # update the best set scores
        @inbounds for varix in used_varixs
            # adjust the set scores by delta
            setix = cur_cover.var2set[varix]
            selix = set2sel(cover_coll, setix)
            @assert selix > 0
            new_score = setscore(cover_coll, setix, cur_cover, selix, varix)
            cur_score = setscore(cover_coll, setix, selix)
            if isfinite(new_score) && (!isfinite(cur_score) || cur_score > new_score)
                scores_updated = true
                cover_coll.sel2cover[selix] = -1 # mark for setting to the cover_pos
            end
        end
        if !scores_updated
            verbose && @info("No global set scores improvement")
            break
        end
        # save the current cover
        insert!(cover_coll.results, cover_pos, cur_cover)
        verbose && @info("Cover collected")
        # update pointers to the best covers for the sets
        for i in eachindex(cover_coll.sel2cover)
            if cover_coll.sel2cover[i] == -1
                cover_coll.sel2cover[i] = cover_pos
            elseif cover_coll.sel2cover[i] >= cover_pos
                # the variant has moved down
                cover_coll.sel2cover[i] += 1
            end
        end
        if params.max_covers > 0 && length(cover_coll) >= params.max_covers
            verbose && @info("Maximal number of covers collected")
            break
        end
        if all(x -> x > 0, cover_coll.sel2cover)
            verbose && @info("All sets assigned to covers")
            break
        end
        # exclude current solution
        cover_problem = exclude_vars(cover_problem,
                            selectvars(cover_problem, mosaic, cur_cover.weights))
    end
    verbose && @info("$(length(cover_coll)) cover(s) collected")
    return cover_coll
end


# DataFrame report for each mask X each set,
# where all set that are covered in at least one mask are considered
# only the best covers used
"""
Convert `covers`, a collection of the covers of `mosaic`, into a `DataFrame`.
"""
function DataFrames.DataFrame(covers::CoverCollection, mosaic::SetMosaic;
                              min_nmasked::Integer=0, min_weight::Real=0.0, best_only::Bool=true)
    # collect all sets that are covered in at least one mask
    selsets = Set{Int}(covers.setixs)
    nmasks = size(covers.elmasks, 2)
    selsize_v = [setsize(mosaic, i) for i in covers.setixs]
    set2sel = Dict(zip(covers.setixs, eachindex(covers.setixs)))
    nmasked_mtx = nmasked_perset(mosaic, covers.elmasks, set2sel)
    coverix_v = Vector{Int}()
    setix_v = Vector{Int}()
    maskix_v = Vector{Int}()
    nmasked_v = Vector{Int}()
    nunmasked_v = Vector{Int}()
    weight_v = Vector{Float64}()
    cover_score_v = Vector{Float64}()
    set_score_enriched_v = Vector{Float64}()
    set_score_overlap_v = Vector{Float64}()

    @inbounds for (coverix, cover) in enumerate(covers.results)
        for (varix, setix) in enumerate(cover.var2set)
            selix = set2sel[setix]
            set_weight = cover.weights[varix]
            set_weight <= min_weight && continue
            set_score = setscore(covers, setix, cover, selix, varix)
            best_only && covers.sel2cover[selix] != coverix && continue
            for maskix in 1:nmasks
                min_nmasked > 0 && nmasked_mtx[selix, maskix] < min_nmasked && continue
                push!(coverix_v, coverix)
                push!(setix_v, setix)
                push!(maskix_v, maskix)
                push!(weight_v, set_weight)
                push!(cover_score_v, cover.agg_total_score)
                push!(set_score_enriched_v, set_score)
                push!(nmasked_v, nmasked_mtx[selix, maskix])
                push!(nunmasked_v, selsize_v[selix] - last(nmasked_v))
                push!(set_score_overlap_v, overlap_score(last(nmasked_v), selsize_v[selix],
                                                   covers.total_masked[maskix], nelements(mosaic),
                                                   1.0 #= ignore relevance =#, covers.cover_params))
            end
        end
    end
    DataFrame(cover_ix = coverix_v,
              set_ix = setix_v,
              set_id = mosaic.ix2set[setix_v],
              mask_ix = maskix_v,
              cover_score = cover_score_v,
              nmasked = nmasked_v,
              nunmasked = nunmasked_v,
              set_relevance = mosaic.set_relevances[setix_v],
              set_weight = weight_v,
              set_score_enriched = set_score_enriched_v,
              set_score_overlap = set_score_overlap_v)
end
