"""
Parameters for the `CoverProblem` (Optimal Enriched-Set Cover).
"""
struct CoverParams
    sel_prob::Float64           # prior probability to select the set, penalizes non-zero weights
    min_weight::Float64         # minimal non-zero set probability
    set_relevance_shape::Float64# how much set relevance affects set score, 0 = no effect
    setXset_factor::Float64     # how much set-set overlaps are penalized (setXset_score scale), 0 = no penalty
    setXset_shape::Float64      # how much set-set overlaps are penalized (setXset_score shape), 0 = no penalty
    maskXmask_factor::Float64   # how much activating overlapping sets in different masks is penalized
    maskXmask_shape::Float64    # how much activating overlapping sets in different masks is penalized

    function CoverParams(; sel_prob::Number=0.5, min_weight::Number = 1E-2,
                         set_relevance_shape::Number=1.0,
                         setXset_factor::Number=1.0, setXset_shape::Number=1.0,
                         maskXmask_factor::Number=1.0, maskXmask_shape::Number=1.0)
        (0.0 < sel_prob <= 1.0) || throw(ArgumentError("`set_prob` must be within (0,1] range"))
        (0.0 < min_weight <= 1.0) || throw(ArgumentError("`min_weight` must be within (0,1] range"))
        (0.0 <= set_relevance_shape) || throw(ArgumentError("`set_relevance_shape` must be ≥0"))
        (0.0 <= setXset_factor) || throw(ArgumentError("`setXset_factor` must be ≥0"))
        (0.0 <= setXset_shape) || throw(ArgumentError("`setXset_shape` must be ≥0"))
        (0.0 <= maskXmask_factor) || throw(ArgumentError("`maskXmask_factor` must be ≥0"))
        (0.0 <= maskXmask_shape) || throw(ArgumentError("`maskXmask_shape` must be ≥0"))
        new(sel_prob, min_weight, set_relevance_shape,
            setXset_factor, setXset_shape,
            maskXmask_factor, maskXmask_shape)
    end
end

"""
Linear component of an individual set score for the `CoverProblem`.
Doesn't take into account the overlap with the other selected sets.
"""
function standalonesetscore(masked::Number, set::Number, total_masked::Number, total::Number, relevance::Number, params::CoverParams)
    # FIXME is it just tail=:both for one set of parameters
    #= P-value for masked-vs-set overlap enriched =# res = logpvalue(masked, set, total_masked, total)*(relevance^params.set_relevance_shape) #-
    #= P-value for unmasked-vs-set overlap enriched =# #logpvalue(set - masked, set, total - total_masked, total)
    @assert !isnan(res) "masked=$masked set=$set total_masked=$total_masked total=$total relevance=$relevance res=NaN"
    return res
end

standalonesetscore(set::MaskedSet, mosaic::MaskedSetMosaic, params::CoverParams) =
    standalonesetscore(set.nmasked, set.nmasked + set.nunmasked,
                       nmasked(mosaic, set.mask), nelements(mosaic),
                       mosaic.original.set_relevances[set.set], params)

"""
Optimal Enriched-Set Cover problem -- choose the sets from the collection to cover
the masked(selected) elements.
The optimal sets cover `C` needs to deliver to goals:
* minimize the P-values of masked elements enrichment for each of `C` sets
* minimize the P-values of the pairwise non-overlap of `C` sets with each other.

Fuzzy set selection is possible -- each set is assigned a weight from `[0, 1]` range.
"""
struct CoverProblem
    params::CoverParams

    var_scores::Vector{Float64}
    varXvar_scores::Matrix{Float64}

    function CoverProblem(params::CoverParams,
                          var_scores::Vector{Float64},
                          varXvar_scores::Matrix{Float64})
        length(var_scores) == size(varXvar_scores, 1) == size(varXvar_scores, 2) ||
            throw(ArgumentError("var_scores and varXvar_scores sizes do not match"))
        new(params, var_scores, varXvar_scores)
    end
end

function CoverProblem(mosaic::MaskedSetMosaic, params::CoverParams = CoverParams())
    var_scores = standalonesetscore.(mosaic.maskedsets, mosaic, params) .- log(params.sel_prob)
    # prepare varXvar scores
    varXvar_scores = zeros(eltype(mosaic.original.setXset_scores), length(var_scores), length(var_scores))
    min_sXs = Inf # minimum finite varXvar_scores element
    k_sXs = params.setXset_factor
    k_mXm = params.maskXmask_factor * params.setXset_factor
    α_sXs = params.setXset_shape
    α_mXm = params.maskXmask_shape * params.setXset_shape
    @inbounds for (i, iset) in enumerate(mosaic.maskedsets)
        for (j, jset) in enumerate(mosaic.maskedsets)
            sXs = mosaic.original.setXset_scores[iset.set, jset.set]
            if iset.set == jset.set
                # no penalty for the overlap with itself, even in different masks
                # (in which case it's zero to encourage the reuse of the same sets to cover different masks)
                vXv = zero(eltype(varXvar_scores))
            elseif !isfinite(sXs)
                vXv = -Inf # fix later with min_score
            else
                if sXs < min_sXs
                    min_sXs = sXs
                end
                if iset.mask == jset.mask
                    vXv = -((1.0-sXs)^α_sXs-1.0)*k_sXs
                else
                    # scale the original setXset score by maskXmask_factor if
                    # the sets are from different masks,
                    vXv = -((1.0-sXs)^α_mXm-1.0)*k_mXm
                end
            end
            varXvar_scores[i, j] = vXv
        end
    end
    # replace infinite varXvar score with the minimal finite setXset score
    sXs_min = -((1.0 - 1.25*min_sXs)^α_sXs-1.0)*k_sXs
    mXm_min = -((1.0 - 1.25*min_sXs)^α_mXm-1.0)*k_mXm
    @inbounds for i in eachindex(varXvar_scores)
        if !isfinite(varXvar_scores[i])
            v1, v2 = ind2sub(size(varXvar_scores), i)
            warn("var[$v1]×var[$v2] score is $(varXvar_scores[i])")
            if varXvar_scores[i] < 0.0
                set1 = mosaic.maskedsets[v1]
                set2 = mosaic.maskedsets[v2]
                varXvar_scores[i] = set1.mask == set2.mask ? sXs_min : mXm_min
            end
        end
    end
    CoverProblem(params, var_scores, varXvar_scores)
end

"""
Total number of problem variables.
"""
nvars(problem::CoverProblem) = length(problem.var_scores)

"""
Construct JuMP quadratic minimization model with linear contraints for the given OESC problem.
"""
function opt_model(problem::CoverProblem)
    m = JuMP.Model()
    nw = nvars(problem)
    @variable(m, 0.0 <= w[1:nw] <= 1.0)
    @objective(m, :Min, dot(vec(problem.var_scores), w) - dot(w, problem.varXvar_scores * w))
    return m
end

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

"""
Score (probability) of the OESC coverage.

* `w` probabilities of the sets being covered
"""
function score(problem::CoverProblem, w::Vector{Float64})
    #pen = fix_cover_weights!(w) # FIXME throw an error?
    dot(problem.var_scores - problem.varXvar_scores * w, w)
end

"""
Result of `optimize(CoverProblem)`.
"""
struct CoverProblemResult
    weights::Vector{Float64}
    var_scores::Vector{Float64}
    total_score::Float64
end

"""
Optimize the cover problem.
"""
function optimize(problem::CoverProblem;
                  ini_weights::Vector{Float64} = rand(nvars(problem)),
                  #iterations::Int = 100,
                  solver::MathProgBase.SolverInterface.AbstractMathProgSolver = default_solver())
    (nvars(problem) == 0) && return CoverProblemResult(Vector{Float64}(), Vector{Float64}(), 0.0)

    # Perform the optimization
    #try
    # using JuMP
    m = opt_model(problem)
    setsolver(m, solver)
    solve(m)
    w = copy(getvalue(getindex(m, :w)))
    # remove small non-zero probabilities due to optimization method errors
    for i in eachindex(w)
        @inbounds if w[i] < problem.params.min_weight
            w[i] = 0.0
        end
    end
    return CoverProblemResult(w, problem.var_scores .* w, getobjectivevalue(m))
    #catch x
    #    warn("Exception in optimize(CoverProblem): $x")
    #    return nothing
    #end
end
