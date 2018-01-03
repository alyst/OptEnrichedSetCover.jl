"""
Parameters for the `CoverProblem` (Optimal Enriched-Set Cover).
"""
struct CoverParams
    overlap_penalty::Float64    # ≥0, how much set overlaps are penalized
    sel_prob::Float64           # prior probability to select the set, penalizes non-zero weights
    min_weight::Float64         # minimal non-zero set probability

    function CoverParams(; overlap_penalty::Number = 1.0, sel_prob::Number = 0.9, min_weight::Number = 1E-2)
        (0.0 <= overlap_penalty) || throw(ArgumentError("`overlap_penalty` must be ≥0"))
        (0.0 < sel_prob <= 1.0) || throw(ArgumentError("`set_prob` must be within (0,1] range"))
        (0.0 < min_weight <= 1.0) || throw(ArgumentError("`min_weight` must be within (0,1) range"))
        new(overlap_penalty, sel_prob, min_weight)
    end
end

"""
Linear component of an individual set score for the `CoverProblem`.
Doesn't take into account the overlap with the other selected sets.
"""
function independentsetscore(masked::Number, set::Number, total_masked::Number, total::Number, params::CoverParams)
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

    set_scores::Matrix{Float64}
    setXset_scores::Matrix{Float64}

    function CoverProblem(mosaic::MaskedSetMosaic, params::CoverParams = CoverParams())
		@inbounds set_scores = Float64[independentsetscore(mosaic.nmasked_perset[i, j], setsize(mosaic, i),
                                                           nmasked(mosaic, j), nelements(mosaic), params) for i in 1:nsets(mosaic), j in 1:nmasks(mosaic)]
        # preprocess setXset scores matrix for numerical solution
        @inbounds setXset_scores = scale!(mosaic.original.setXset_scores[mosaic.setixs, mosaic.setixs],
                                          params.overlap_penalty)
        # replace infinite setXset score with the minimal finite setXset score
        min_score = 0.0
        @inbounds for i in eachindex(setXset_scores)
            if !isfinite(setXset_scores[i])
                s1, s2 = ind2sub(size(setXset_scores), i)
                warn("set[$s1]×set[$s2] score is $(setXset_scores[i])")
            elseif setXset_scores[i] < min_score
                min_score = setXset_scores[i]
            end
        end
        @inbounds for i in eachindex(setXset_scores)
            if isinf(setXset_scores[i]) && setXset_scores[i] < 0.0
                setXset_scores[i] = min_score
            end
        end
        # activating the same set in another mask doesn't introduce the penalty
        @inbounds for i in 1:size(setXset_scores, 1)
            setXset_scores[i, i] = 0
        end
        # activating the set in a given mask introduces the penalty
        # (too constrain the number of activated sets)
        const log_selp = log(params.sel_prob)
        multi_setXset_scores = repeat(setXset_scores, outer=[nmasks(mosaic), nmasks(mosaic)])
        @inbounds for i in 1:size(multi_setXset_scores, 1)
            multi_setXset_scores[i, i] = log_selp
        end
        new(params, set_scores, multi_setXset_scores)
    end
end

"""
Total number of sets in the collection.
"""
nsets(problem::CoverProblem) = size(problem.set_scores, 1)

"""
Total number of masks in the collection.
"""
nmasks(problem::CoverProblem) = size(problem.set_scores, 2)

nweights(problem::CoverProblem) = length(problem.set_scores)

"""
Construct JuMP quadratic minimization model with linear contraints for the given OESC problem.
"""
function opt_model(problem::CoverProblem)
    m = JuMP.Model()
    nw = nweights(problem)
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
    weights::Matrix{Float64}
    score::Float64

    CoverProblemResult(weights::Matrix{Float64}, score::Float64) =
        new(weights, score)
end

"""
Optimize the cover problem.
"""
function optimize(problem::CoverProblem;
                  ini_weights::Vector{Float64} = rand(nweights(problem)),
                  #iterations::Int = 100,
                  solver::MathProgBase.SolverInterface.AbstractMathProgSolver = IpoptSolver(print_level=0))
    (nweights(problem) == 0) && return CoverProblemResult(Matrix{Float64}(0,0), 0.0)

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
    return CoverProblemResult(reshape(w, (nsets(problem), nmasks(problem))), getobjectivevalue(m))
    #catch x
    #    warn("Exception in optimize(CoverProblem): $x")
    #    return nothing
    #end
end
