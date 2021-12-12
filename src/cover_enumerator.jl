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
struct CoverCollection
    cover_params::CoverParams
    enum_params::CoverEnumerationParams
    total_masked::Vector{Int}         # total masked elements in the mosaic
    ix2experiment::Vector             # vector of experiment Ids, copied from masked mosaic
    elmasks::BitMatrix                # FIXME the elements mask, a workaround to avoid copying the whole mosaic upon serialization
    setixs::Vector{Int}               # sets referenced by the MaskedSetMosaic
    base_selscores::Vector{Float64}   # base selected set scores
    sel2cover::Vector{Int}            # best-scoring cover for the given selected set
    results::Vector{MultiobjCoverProblemResult}

    function CoverCollection(problem::MultiobjCoverProblem, mosaic::MaskedSetMosaic,
                             params::CoverEnumerationParams)
        # check the problem is compatible with the mosaic
        nvars(problem) == nsets(mosaic) ||
                throw(ArgumentError("CoverProblem is not compatible to the MaskedSetMosaic: number of vars and sets differ"))
        #nmasks(problem) == nmasks(mosaic) || throw(ArgumentError("CoverProblem is not compatible to the MaskedSetMosaic: number of masks differ"))
        new(problem.params, params, mosaic.total_masked, mosaic.ix2experiment, mosaic.elmasks,
            problem.var2set, copy(problem.var_scores), zeros(Int, nvars(problem)),
            Vector{MultiobjCoverProblemResult}())
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

function setscore(covers::CoverCollection, setix::Integer,
                  cover::MultiobjCoverProblemResult,
                  selix::Int=0, solix::Int=0, varix::Int=0)
    selix = set2sel(covers, setix, selix)
    (selix == 0) && return NaN # the set is not a part of the cover collection
    if solix == 0
        solix = best_index(cover)
    end
    if varix == 0
        varix = set2var(cover, setix)
        (varix == 0) && return NaN # the set is not a part of the cover problem
    elseif cover.var2set[varix] != setix
        throw(ArgumentError("Incorrect varix hint ($varix)"))
    end
    varw = cover.varweights[varix, solix]
    if varw > 0.0
        # problem set score + delta score for the best variant, where it was covered - log(set weight)
        return covers.base_selscores[selix] - log(varw)
    else
        # the set is not selected
        return Inf
    end
end

function setscore(covers::CoverCollection, setix::Int, selix::Int=0, solix::Int=0)
    selix = set2sel(covers, setix, selix)
    (selix == 0) && return NaN # the set is not a part of the cover collection
    if covers.sel2cover[selix] > 0
        # problem set score + delta score for the best variant, where it was covered - log(set weight)
        return setscore(covers, setix, covers.sel2cover[selix], selix, solix)
    else
        # the set is not covered
        return NaN
    end
end

nexperiments(covers::CoverCollection) = size(covers.elmasks, 2)
Base.length(covers::CoverCollection) = length(covers.results)
Base.isempty(covers::CoverCollection) = isempty(covers.results)

function OptimizerParams(problem_type::Symbol, params::ParamsDict)
    if problem_type == :quadratic
        return QuadraticOptimizerParams(;params...)
    elseif problem_type == :multiobjective
        return MultiobjOptimizerParams(;params...)
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
    append!(CoverCollection(cover_problem, mosaic, enum_params),
            cover_problem, mosaic, enum_params, opt_params, verbose)
end

# append the optimal covers of `cover_problem` into `cover_coll`
function Base.append!(cover_coll::CoverCollection,
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
        if isempty(cur_cover)
            verbose && @info("No covers")
            break
        end
        verbose && @info("New cover found (score=$(cur_cover.scores[best_index(cur_cover)]), aggscore=$(best_aggscore(cur_cover))), processing...")
        used_varixs = findall(w -> w > cover_problem.params.min_weight, best_varweights(cur_cover))
        if isempty(used_varixs)
            verbose && @info("Cover is empty")
            break
        end
        # store cover in the list maintaining sorting order
        cover_pos = searchsortedlast(cover_coll.results, cur_cover, by=best_aggscore)+1
        if cover_pos > 1
            cover = cover_coll.results[cover_pos-1]
            if abs(best_aggscore(cur_cover) - best_aggscore(cover)) <= score_threshold &&
               cur_cover.var2set == cover.var2set #= should never happen =#
                w_cur = best_varweights(cur_cover)
                w_cover = best_varweights(cover)
                if all(i -> (@inbounds return abs(w_cur[i] - w_cover[i]) <= weight_threshold),
                    eachindex(w_cur, w_cover))
                    verbose && @info("Duplicate solution")
                    break
                end
            end
        end
        if cover_pos > 1 && isa(params.max_cover_score_delta, Float64)
            # not the best cover
            cover_score_delta = best_aggscore(cur_cover) - best_aggscore(cover_coll.results[1])
            if cover_score_delta > params.max_cover_score_delta
                verbose && @info("Cover score_delta=$(cover_score_delta) above threshold")
                break
            end
        end
        if isa(params.max_set_score, Float64)
            max_set_score = params.max_set_score::Float64
            if all(i -> (@inbounds return all(x -> x > max_set_score, view(cover_problem.var_scores, i, :))),
                    used_varixs)
                verbose && @info("All set scores above $(max_set_score)")
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
            new_score = setscore(cover_coll, setix, cur_cover, selix, 0, varix)
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
                            selectvars(cover_problem, mosaic, best_varweights(cur_cover)))
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
                              min_nmasked::Integer=0, min_weight::Real=0.0,
                              best_only::Bool=true, params::Union{CoverParams, Nothing} = nothing)
    # collect all sets that are covered in at least one mask
    selsets = Set{Int}(covers.setixs)
    nexps = nexperiments(covers)
    selsize_v = setsize.(Ref(mosaic), covers.setixs)
    set2sel = Dict(zip(covers.setixs, eachindex(covers.setixs)))
    nmasked_mtx = nmasked_perset(mosaic, covers.elmasks, set2sel)
    coverix_v = Vector{Int}()
    setix_v = Vector{Int}()
    expix_v = Vector{Int}()
    nmasked_v = Vector{Int}()
    nunmasked_v = Vector{Int}()
    weight_v = Vector{Float64}()
    cover_total_score_v = Vector{Float64}()
    set_cover_score_v = Vector{Float64}()
    set_overlap_logpvalue_v = Vector{Float64}()

    @inbounds for (coverix, cover) in enumerate(covers.results)
        sol_ix = best_index(cover, params)
        sol_aggscore = params === nothing ?
                       best_aggscore(cover) :
                       aggscore(cover.scores[sol_ix], params)
        sol_varweights = varweights(cover, sol_ix)
        for (varix, setix) in enumerate(cover.var2set)
            selix = set2sel[setix]
            set_weight = sol_varweights[varix]
            set_weight <= min_weight && continue
            set_score = setscore(covers, setix, cover, selix, sol_ix, varix)
            best_only && covers.sel2cover[selix] != coverix && continue
            for expix in 1:nexps
                min_nmasked > 0 && nmasked_mtx[selix, expix] < min_nmasked && continue
                push!(coverix_v, coverix)
                push!(setix_v, setix)
                push!(expix_v, expix)
                push!(weight_v, set_weight)
                push!(cover_total_score_v, sol_aggscore)
                push!(set_cover_score_v, set_score)
                push!(nmasked_v, nmasked_mtx[selix, expix])
                push!(nunmasked_v, selsize_v[selix] - last(nmasked_v))
                push!(set_overlap_logpvalue_v, logpvalue(last(nmasked_v), selsize_v[selix],
                                                         covers.total_masked[expix], nelements(mosaic)))
            end
        end
    end
    DataFrame(cover_ix = coverix_v,
              set_ix = setix_v,
              set_id = mosaic.ix2set[setix_v],
              experiment_ix = expix_v,
              experiment_id = covers.ix2experiment[expix_v],
              cover_total_score = cover_total_score_v,
              nmasked = nmasked_v,
              nunmasked = nunmasked_v,
              set_relevance = mosaic.set_relevances[setix_v],
              set_weight = weight_v,
              set_cover_score = set_cover_score_v,
              set_overlap_log10pvalue = set_overlap_logpvalue_v ./ log(10))
end
