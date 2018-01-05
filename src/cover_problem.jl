"""
Parameters for the `CoverProblem` (Optimal Enriched-Set Cover).
"""
struct CoverParams
    sel_prob::Float64           # prior probability to select the set, penalizes non-zero weights
    min_weight::Float64         # minimal non-zero set probability
    setXset_factor::Float64     # how much set-set overlaps are penalized (setXset_score scale), 0 = no penalty
    maskXmask_factor::Float64   # how much activating overlapping sets in different masks is penalized

    function CoverParams(; sel_prob::Number=0.9, min_weight::Number = 1E-2,
                         setXset_factor::Number=1.0, maskXmask_factor::Number=1.0)
        (0.0 < sel_prob <= 1.0) || throw(ArgumentError("`set_prob` must be within (0,1] range"))
        (0.0 < min_weight <= 1.0) || throw(ArgumentError("`min_weight` must be within (0,1] range"))
        (0.0 <= setXset_factor) || throw(ArgumentError("`setXset_factor` must be ≥0"))
        (0.0 <= maskXmask_factor <= 1.0) || throw(ArgumentError("`maskXmask_factor` must be within [0,1] range"))
        new(sel_prob, min_weight, setXset_factor, maskXmask_factor)
    end
end

"""
Linear component of an individual set score for the `CoverProblem`.
Doesn't take into account the overlap with the other selected sets.
"""
function standalonesetscore(masked::Number, set::Number, total_masked::Number, total::Number, params::CoverParams)
    # FIXME is it just tail=:both for one set of parameters
    #= P-value for masked-vs-set overlap enriched =# res = logpvalue(masked, set, total_masked, total) #-
    #= P-value for unmasked-vs-set overlap enriched =# #logpvalue(set - masked, set, total - total_masked, total)
    @assert !isnan(res) "masked=$masked set=$set total_masked=$total_masked total=$total res=NaN"
    return res
end

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

    set_scores::Vector{Float64}
    setXset_scores::Matrix{Float64}

    function CoverProblem(params::CoverParams,
                          set_scores::Vector{Float64},
                          setXset_scores::Matrix{Float64})
        length(set_scores) == size(setXset_scores, 1) == size(setXset_scores, 2) ||
            throw(ArgumentError("set_scores and setXset_scores sizes do not match"))
        new(params, set_scores, setXset_scores)
    end
end

function CoverProblem(mosaic::MaskedSetMosaic, params::CoverParams = CoverParams())
    # activating the set in a given mask introduces the penalty
    # (to constrain the number of activated sets)
    const log_selp = log(params.sel_prob)
	@inbounds set_scores = Float64[standalonesetscore(s.nmasked, s.nmasked + s.nunmasked,
                                                      nmasked(mosaic, s.mask), nelements(mosaic), params) - log_selp for s in mosaic.maskedsets]
    # prepare setXset scores
    setXset_scores = zeros(eltype(mosaic.original.setXset_scores), length(set_scores), length(set_scores))
    min_score = Inf # minimum finite setXset_scores element
    @inbounds for (i, iset) in enumerate(mosaic.maskedsets)
        for (j, jset) in enumerate(mosaic.maskedsets)
            sc = params.setXset_factor * mosaic.original.setXset_scores[iset.set, jset.set]
            if iset.mask != jset.mask
                if iset.set != jset.set
                    # scale the original setXset score by maskXmask_factor if
                    # the sets are from different masks,
                    sc *= params.maskXmask_factor
                else
                    # unless (i, j) point to the same set
                    # (in which case it's zero to encourage the reuse of the same sets to cover different masks)
                    sc = zero(eltype(setXset_scores))
                end
            end
            setXset_scores[i, j] = sc
            if isfinite(sc) && sc < min_score
                min_score = sc
            end
        end
    end
    # replace infinite setXset score with the minimal finite setXset score
    @inbounds for i in eachindex(setXset_scores)
        if !isfinite(setXset_scores[i])
            s1, s2 = ind2sub(size(setXset_scores), i)
            warn("set[$s1]×set[$s2] score is $(setXset_scores[i])")
            if setXset_scores[i] < 0.0
                setXset_scores[i] = min_score
            end
        end
    end
    CoverProblem(params, set_scores, setXset_scores)
end

"""
Total number of problem variables.
"""
nvars(problem::CoverProblem) = length(problem.set_scores)

"""
Construct JuMP quadratic minimization model with linear contraints for the given OESC problem.
"""
function opt_model(problem::CoverProblem)
    m = JuMP.Model()
    nw = nvars(problem)
    @variable(m, 0.0 <= w[1:nw] <= 1.0)
    @objective(m, :Min, dot(vec(problem.set_scores), w) - dot(w, problem.setXset_scores * w))
    return m
end

function fix_uncov_probs!(uncov_probs::Matrix{Float64})
    pen = 0.0
    @inbounds for i in eachindex(uncov_probs)
        prob = uncov_probs[i]
        prob_new = clamp(prob, 0.0, 1.0)
        pen += abs2(prob_new - prob)
        uncov_probs[i] = prob_new
    end
    return pen
end

"""
Score (probability) of the OESC coverage.

* `w` probabilities of the sets being covered
"""
function score(problem::CoverProblem, w::Vector{Float64})
    # FIXME throw an error?
    #pen = fix_uncov_probs!(uncov_probs)
    dot(problem.set_scores - problem.setXset_scores * w, w)
end

"""
Result of `optimize(CoverProblem)`.
"""
struct CoverProblemResult
    weights::Vector{Float64}
    total_score::Float64
end

"""
Optimize the cover problem.
"""
function optimize(problem::CoverProblem;
                  ini_weights::Vector{Float64} = rand(nvars(problem)),
                  #iterations::Int = 100,
                  solver::MathProgBase.SolverInterface.AbstractMathProgSolver = IpoptSolver(print_level=0))
    (nvars(problem) == 0) && return CoverProblemResult(Vector{Float64}(), 0.0)

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
    return CoverProblemResult(w, getobjectivevalue(m))
    #catch x
    #    warn("Exception in optimize(CoverProblem): $x")
    #    return nothing
    #end
end
