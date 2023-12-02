"""
    CoverParams(; [sel_tax=0.0], [set_shape=1.0],
                  [min_weight=1E-2], [mask_discount=0.9], [setXset_factor=1.0],
                  [uncovered_factor=0.1], [covered_factor=0.001],
                  [set_relevance_shape=0.5], [set_relevance_min=0.5]) -> CoverParams

Specify parameters for the [`AbstractCoverProblem`](@ref).
One can specify a non-default value for a particular parameter.

See [cover quality score](@ref cover_quality) for the detailed description of the parameters.

# Arguments
  * `sel_tax`: ``\\epsilon``, the constant added to the set score of each selected set to penalize
     the number of sets in the cover
  * `set_shape`: ``\\beta``, applied to set or set×set scores
  * `min_weight`: minimal non-zero set probability ???
  * `mask_discount`: ``\\alpha``, how much the overlap score of each subsequent mask (from most to less enriched)
    is discounted
  * `setXset_factor`: ``w_r``, the weight of *redundancy* [cover quality](@ref cover_quality) component
    (`setXset_score` scale), 0 = no redunancy penalty
  * `uncovered_factor`: ``w_u``, the weight of *uncovered hits* component in [cover quality](@ref cover_quality)
  * `covered_factor`: ``w_c``, the weight of *covered non-hits* component in [cover quality](@ref cover_quality)
  * `set_relevance_shape`: ``\\beta_L``, how much set relevance affects set score, 0 = no effect
  * `set_relevance_min`: ``L_{\\min``, if shaped relevance is below, it's set to `set_relevance_min`
"""
struct CoverParams
    sel_tax::Float64
    set_shape::Float64
    min_weight::Float64
    mask_discount::Float64
    setXset_factor::Float64
    uncovered_factor::Float64
    covered_factor::Float64
    set_relevance_shape::Float64
    set_relevance_min::Float64

    function CoverParams(;
                         sel_tax::Real=0.0, set_shape::Real=1.0,
                         min_weight::Real=1E-2, mask_discount::Real=0.9,
                         setXset_factor::Real=1.0,
                         uncovered_factor::Real=0.1, covered_factor::Real=0.001,
                         set_relevance_shape::Real=0.5,
                         set_relevance_min::Real=0.5)
        (0.0 < min_weight <= 1.0) || throw(ArgumentError("`min_weight` must be within (0,1] range"))
        (0.0 <= mask_discount <= 1.0) || throw(ArgumentError("`mask_discount` must be within [0,1] range"))
        (0.0 <= set_shape) || throw(ArgumentError("`set_shape` must be ≥0"))
        (0.0 <= setXset_factor) || throw(ArgumentError("`setXset_factor` must be ≥0"))
        (0.0 <= uncovered_factor) || throw(ArgumentError("`uncovered_factor` must be ≥0"))
        (0.0 <= covered_factor) || throw(ArgumentError("`covered_factor` must be ≥0"))
        (0.0 <= set_relevance_shape) || throw(ArgumentError("`set_relevance_shape` must be ≥0"))
        (0.0 <= set_relevance_min <= 1) || throw(ArgumentError("`set_relevance_min` must be within [0, 1] range"))
        new(sel_tax, set_shape, min_weight, mask_discount,
            setXset_factor, uncovered_factor, covered_factor,
            set_relevance_shape, set_relevance_min)
    end
end

"""
Linear component of an individual masked set score for the `CoverProblem`.
Doesn't take into account the overlap with the other selected sets.
"""
function overlap_score(masked::Number, set::Number, total_masked::Number, total::Number,
                       relevance::Number, params::CoverParams)
    # FIXME is it just tail=:both for one set of parameters
    #= P-value for masked-vs-set overlap enriched =# res = -(-logpvalue(masked, set, total_masked, total))^params.set_shape*max(relevance^params.set_relevance_shape, params.set_relevance_min) #-
    #= P-value for unmasked-vs-set overlap enriched =# #logpvalue(set - masked, set, total - total_masked, total)
    @assert !isnan(res) "masked=$masked set=$set total_masked=$total_masked total=$total relevance=$relevance res=NaN"
    return res
end

function _setweight(mosaic::AbstractWeightedSetMosaic,
                    setix::Integer, expix::Integer,
                    relevance::Number, params::CoverParams;
                    filter::Bool = true)
    w = _setweight(mosaic, setix, expix, filter=filter)
    ismissing(w) && return w
    res = -(-w)^params.set_shape *
            max(relevance^params.set_relevance_shape, params.set_relevance_min) #-
    #= P-value for unmasked-vs-set overlap enriched =# #logpvalue(set - masked, set, total - total_masked, total)
    @assert !isnan(res) "masked=$masked set=$set total_masked=$total_masked total=$total relevance=$relevance res=NaN"
    return res
end

var2set(mosaic::AbstractWeightedSetMosaic) = eachindex(mosaic.loc2glob_setix)

# calculates var scores
# var score is:
# the sum of overlap scores with all masked sets
# minus the log(var selection probability)
# plus nets*log(nsets), where nsets arr all masked sets overlapping with var
function var_scores(mosaic::AbstractWeightedSetMosaic, params::CoverParams,
                    var2setix::AbstractVector{Int} = var2set(mosaic))
    # find relevant tiles and
    # sum masked elements per tile in all masks and unmasked elements per tile in union mask
    # note: in masks there could be element not covered by any set/tile, so these counts refer to coverable elements
    # calculate the sum of scores of given set in each mask
    v_scores = Vector{Float64}(undef, length(var2setix))
    scores = Vector{Float64}() # temp array containing var overlap scores for each mask it overlaps

    @inbounds for (varix, setix) in enumerate(var2setix)
        empty!(scores)
        glob_setix = mosaic.loc2glob_setix[setix]
        setrelevance = originalmosaic(mosaic).set_relevances[glob_setix]
        for expix in 1:nexperiments(mosaic)
            w = _setweight(mosaic, setix, expix, setrelevance, params)
            ismissing(w) || push!(scores, w)
        end
        sort!(scores)
        scoresum = 0.0
        w = 1.0
        @inbounds for score in scores
            scoresum += score * w
            w *= params.mask_discount
        end
        v_scores[varix] = scoresum + params.sel_tax
    end
    return v_scores
end

# prepares
# 1) tile-to-var mapping (SparseMaskMatrix, only relevant tiles are considered)
# 2) var-to-(tile-in-mask) mapping (SparseMatrixCSC, only relevant tiles are considered)
# 3) the total number of masked elements in all masks for each tile
# 4) the total number of unmasked elements in all active masks (overlapping with the tile) for each tile
function varXtiles(mosaic::MaskedSetMosaic,
                   var2setix::AbstractVector{Int},
                   params::CoverParams)
    # find relevant tiles and
    # sum masked elements per tile in all masks and unmasked elements per tile in union mask
    # note: in masks there could be element not covered by any set/tile, so these counts refer to coverable elements
    tileixs = Vector{Int}()
    tile2nmasked = Vector{Int}()
    tile2nunmasked = Vector{Int}()
    vXmt_els = Vector{Tuple{Int, Int, Int, Float64}}()
    # max rank of the setXexp overlap that is still relevant
    setXexp_olaps = fill(0.0, nexperiments(mosaic))

    @inbounds for (varix, setix) in enumerate(var2setix)
        globsetix = mosaic.loc2glob_setix[setix]
        for expix in 1:nexperiments(mosaic)
            setXexp_olaps[expix] = coalesce(_setweight(mosaic, setix, expix), Inf)
        end
        # rank 1 = most significant (weights are negative)
        setXexp_olap_ranks = denserank(setXexp_olaps)

        set_tiles = view(mosaic.original.tileXset, :, globsetix)
        for tileix in set_tiles
            tile_elms = view(mosaic.original.elmXtile, :, tileix)
            tilepos = searchsortedfirst(tileixs, tileix)
            if tilepos > length(tileixs) || tileixs[tilepos] != tileix # do this once per tile
                insert!(tileixs, tilepos, tileix)
                insert!(tile2nunmasked, tilepos, length(tile_elms) - sum(view(mosaic.elunionmask, tile_elms)))
                nmasked = 0
                for (expix, olap) in enumerate(view(mosaic.setXexp_olaps, setix, :))
                    if !ismissing(olap)
                        nmasked += sum(view(mosaic.elmasks, tile_elms, expix))
                    end
                end
                insert!(tile2nmasked, tilepos, nmasked)
            end
            for (expix, olap) in enumerate(view(mosaic.setXexp_olaps, setix, :))
                ismissing(olap) && continue
                nmaskxtile = sum(view(mosaic.elmasks, tile_elms, expix))
                olap_rank = setXexp_olap_ranks[expix]
                olap_w = params.mask_discount ^ (olap_rank - 1)
                if nmaskxtile < length(tile_elms) && # there are unmasked tile elements
                   olap_w >= 0.1                     # the mask has relevant overlap with the var
                    push!(vXmt_els, (varix, expix, tileix, olap_w * (length(tile_elms) - nmaskxtile)))
                end
            end
        end
    end

    tileXvar = mosaic.original.tileXset[tileixs, mosaic.loc2glob_setix[var2setix]]
    if !isempty(tileXvar)
        # try to optimize tileXvar: find tiles that refer to the identical set of vars and merge their stats
        varXtile = transpose(tileXvar)
        vars2tiles = Vector{Tuple{Vector{Int}, Vector{Int}}}()
        for tileix in axes(varXtile, 2)
            varixs = varXtile[:, tileix]
            v2tpos = searchsortedfirst(vars2tiles, (varixs, 0), by=first)
            if v2tpos > length(vars2tiles) || first(vars2tiles[v2tpos]) != varixs
                insert!(vars2tiles, v2tpos, (varixs, [tileix]))
            else
                push!(vars2tiles[v2tpos][2], tileix)
            end
        end
        if length(vars2tiles) < size(tileXvar, 1) # there are tiles to merge
            # take the 1st (it's also minimal since tileixs are sorted) of each merged tiles
            tileXvar = tileXvar[[first(tileixs) for (_, tileixs) in vars2tiles], :]
            # aggregate n[un]masked over duplicated tiles
            tile2nmasked = [sum(i -> tile2nmasked[i], tileixs) for (_, tileixs) in vars2tiles]
            tile2nunmasked = [sum(i -> tile2nunmasked[i], tileixs) for (_, tileixs) in vars2tiles]
        end
    end
    if true # optimize varXmaskxtile: find maskXtile that refers to the identical set of vars and sum their vals
        maskxtile2vars = Dict{Tuple{Int, Int}, Vector{Int}}()
        for (varix, expix, tileix, val) in vXmt_els
            varixs = get!(() -> Vector{Int}(), maskxtile2vars, (expix, tileix))
            push!(varixs, varix)
        end
        # each varixs should be sorted, because they are traversed in asc order when constructing vXmt_els
        vars2maskxtiles = Vector{Tuple{Vector{Int}, Vector{Tuple{Int, Int}}}}()
        for (maskxtile, varixs) in pairs(maskxtile2vars)
            v2mtpos = searchsortedfirst(vars2maskxtiles, (varixs, 0), by=first)
            if v2mtpos > length(vars2maskxtiles) || first(vars2maskxtiles[v2mtpos]) != varixs
                insert!(vars2maskxtiles, v2mtpos, (varixs, [maskxtile]))
            else
                # FIXME use OrderedDict to avoid searchsorted?
                maskxtiles = last(vars2maskxtiles[v2mtpos])
                insert!(maskxtiles, searchsortedfirst(maskxtiles, maskxtile), maskxtile)
            end
        end
        maskxtile2ix = Dict{Tuple{Int, Int}, Int}()
        for (i, (_, maskxtiles)) in enumerate(vars2maskxtiles)
            for maskxtile in maskxtiles
                maskxtile2ix[maskxtile] = i
            end
        end
        # replace expix by maskxtile index in vxmxt_els
        @inbounds for (i, (varix, expix, tileix, v)) in enumerate(vXmt_els)
            vXmt_els[i] = (varix, maskxtile2ix[(expix, tileix)], 0, v)
        end
        # sort vxmxt_els by varix, the by maskxtile index to match SpraseMatrixCSC order
        sort!(vXmt_els, lt=(a,b) -> (a[1] < b[1]) || ((a[1] == b[1]) && (a[2] < b[2])))
        # construct sparse matrix
        colptrs = Vector{Int}(undef, length(var2setix) + 1)
        rowvals = Vector{Int}(undef, mapreduce(x -> length(first(x)), +, vars2maskxtiles, init=0))
        nzvals = fill!(Vector{Float64}(undef, length(rowvals)), 0)
        last_coord = (0, 0)
        i = 0
        @inbounds for (varix, maskxtileix, _, v) in vXmt_els
            if last_coord != (varix, maskxtileix)
                @assert i < length(rowvals)
                i += 1
                if last_coord[1] != varix
                    @assert last_coord[1] < varix
                    colptrs[last_coord[1]+1:varix] .= i
                else
                    @assert maskxtileix > last_coord[2]
                end
                last_coord = (varix, maskxtileix)
            end
            rowvals[i] = maskxtileix
            nzvals[i] += v
        end
        @assert i == length(rowvals)
        colptrs[last_coord[1]+1:end] .= length(rowvals) + 1
        maskxtileXvar = SparseMatrixCSC{Float64, Int}(length(vars2maskxtiles), length(var2setix),
                                                      colptrs, rowvals, nzvals)
    end
    return tileXvar, maskxtileXvar, tile2nmasked, tile2nunmasked
end

varXvar_score(setXset::Real, params::CoverParams, scale::Bool = false) =
    ifelse(scale, -(-setXset)^params.set_shape * params.setXset_factor, setXset)

function varXvar_scores(mosaic::AbstractWeightedSetMosaic, params::CoverParams,
                        var2setix::AbstractVector{Int} = var2set(mosaic),
                        scale::Bool = false)
    glob_setixs = mosaic.loc2glob_setix[var2setix]
    vXv_scores = varXvar_score.(mosaic.original.setXset_scores[glob_setixs, glob_setixs], Ref(params), scale)
    for i in eachindex(var2setix) # don't consider self-intersections
        @inbounds vXv_scores[i, i] = zero(eltype(vXv_scores))
    end
    # find minimum finite varXvar_scores element
    vXv_min = Inf
    @inbounds for vXv in vXv_scores
        if isfinite(vXv) && vXv < vXv_min
            vXv_min = vXv
        end
    end
    # replace infinite varXvar score with the minimal finite setXset score
    vXv_subst = varXvar_score(1.25*vXv_min, params, true)
    vXv_cartixs = CartesianIndices(vXv_scores)
    @inbounds for i in eachindex(vXv_scores)
        if !isfinite(vXv_scores[i])
            v1, v2 = Tuple(vXv_cartixs[i])
            @warn "var[$v1]×var[$v2] score is $(vXv_scores[i])"
            if vXv_scores[i] < 0.0
                vXv_scores[i] = vXv_subst
            end
        end
    end
    return vXv_scores
end

# detailed cover score:
# 1 => covered set scores,
# 2 => covered setXset penalties,
# 3 => uncovered masked elements penalties
# 4 => covered unmasked elements penalties
const DetailedScore = NTuple{4, Float64}

aggscore(s::DetailedScore, params::CoverParams) =
    s[1] +
    params.setXset_factor * s[2] +
    # ignore NaNs/Infs for uncovered/covered score components if their factors are 0
    ifelse(params.uncovered_factor != 0, params.uncovered_factor * s[3], zero(s[3])) +
    ifelse(params.covered_factor != 0, params.covered_factor * s[4], zero(s[4]))

"""
    AbstractCoverProblem{T}

Optimal Enriched-Set Cover problem -- choose the sets from the collection `𝒞` to cover
the masked(selected) elements `M`.
The optimal sets cover `C = {c₁, c₂, ..., cₙ} ⊂ 𝒞` has to deliver 3 goals:
* be relevant (i.e. minimize the P-values of `M` and `cᵢ` sets overlap)
* be minimal (i.e. minimize the number of sets in `C`)
* be non-redundant (i.e. minimize the P-values of the pairwise non-overlap of `C` sets with each other).

Fuzzy set selection is possible -- each set is assigned a weight from `[0, 1]` range.
"""
abstract type AbstractCoverProblem{T} end

"""
Total number of problem variables.
"""
nvars(problem::AbstractCoverProblem) = length(problem.var_scores)

# clamps weights to 0..1 range and returns sum(|w - fixed_w|)
function fix_cover_weights!(weights::Matrix{Float64})
    pen = 0.0
    @inbounds for i in eachindex(weights)
        w = uncov_probs[i]
        w_new = clamp(w, 0.0, 1.0)
        pen += abs2(w_new - w)
        weights[i] = w_new
    end
    return pen
end

__check_vars(w::AbstractVector{Float64}, problem::AbstractCoverProblem) =
    length(w) == nvars(problem) ||
        throw(DimensionMismatch("Incorrect length of parameters vector: $(length(w)) ($(nvars(problem)) expected)"))

function selectvars(problem::AbstractCoverProblem,
                    _::AbstractWeightedSetMosaic,
                    weights::AbstractVector{Float64})
    @assert length(weights) == nvars(problem)
    return findall(>(problem.params.min_weight), weights) # > to make it work with min_weight=0
end

abstract type AbstractOptimizerParams{P <: AbstractCoverProblem} end;

problemtype(_::AbstractOptimizerParams{P}) where P = P
