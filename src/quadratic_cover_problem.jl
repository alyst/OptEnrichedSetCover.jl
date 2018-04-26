"""
Quadratic Programming Optimal Enriched-Set Cover problem.
"""
struct QuadraticCoverProblem <: AbstractCoverProblem{Float64}
    params::CoverParams

    var_scores::Vector{Float64}
    varXvar_scores::Matrix{Float64}

    function QuadraticCoverProblem(params::CoverParams,
                          var_scores::Vector{Float64},
                          varXvar_scores::Matrix{Float64})
        length(var_scores) == size(varXvar_scores, 1) == size(varXvar_scores, 2) ||
            throw(ArgumentError("var_scores and varXvar_scores sizes do not match"))
        new(params, var_scores, varXvar_scores)
    end
end

struct QuadraticOptimizerParams <: AbstractOptimizerParams{QuadraticCoverProblem}
    solver::MathProgBase.SolverInterface.AbstractMathProgSolver

    QuadraticOptimizerParams(;
        solver::MathProgBase.SolverInterface.AbstractMathProgSolver = default_quadratic_solver(),
        #iterations::Int = 100,
        #ini_weights::Vector{Float64}
    ) = new(solver)
end

quadratic_solver(opt_params::QuadraticOptimizerParams) = opt_params.solver

function QuadraticCoverProblem(mosaic::MaskedSetMosaic, params::CoverParams = CoverParams())
    var_scores = msetscore_detached.(mosaic.maskedsets, mosaic, params) .- log(params.sel_prob)
    # prepare varXvar scores
    varXvar_scores = zeros(eltype(mosaic.original.setXset_scores), length(var_scores), length(var_scores))
    min_sXs = Inf # minimum finite varXvar_scores element
    @inbounds for (i, iset) in enumerate(mosaic.maskedsets)
        for (j, jset) in enumerate(mosaic.maskedsets)
            sXs = mosaic.original.setXset_scores[iset.set, jset.set]
            if iset.set != jset.set && isfinite(sXs) && sXs < min_sXs
                min_sXs = sXs
            end
            varXvar_scores[i, j] = varXvar_score(sXs, iset, jset, params, true)
        end
    end
    # replace infinite varXvar score with the minimal finite setXset score
    sXs_min = varXvar_score(1.25*min_sXs, MaskedSet(1, 1, 0, 0), MaskedSet(2, 1, 0, 0), params, true)
    mXm_min = varXvar_score(1.25*min_sXs, MaskedSet(1, 1, 0, 0), MaskedSet(2, 2, 0, 0), params, true)
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
    QuadraticCoverProblem(params, var_scores, varXvar_scores)
end

"""
Construct JuMP quadratic minimization model with linear contraints for the given OESC problem.
"""
function opt_model(problem::QuadraticCoverProblem)
    m = JuMP.Model()
    nw = nvars(problem)
    @variable(m, 0.0 <= w[1:nw] <= 1.0)
    @objective(m, :Min, dot(vec(problem.var_scores), w) - dot(w, problem.varXvar_scores * w))
    return m
end

"""
Score of the OESC coverage.

* `w` probabilities of the sets being covered
"""
function score(problem::QuadraticCoverProblem, w::Vector{Float64})
    #pen = fix_cover_weights!(w) # FIXME throw an error?
    dot(problem.var_scores - problem.varXvar_scores * w, w)
end

aggscore(problem::QuadraticCoverProblem, w::Vector{Float64}) = score(w, problem)

"""
    optimize(problem::AbstractCoverProblem, method_params::QuadraticMethodParams)
    optimize(problem::AbstractCoverProblem, [method_kwparams...])

Optimize the cover problem.
"""
function optimize(problem::QuadraticCoverProblem,
                  opt_params::QuadraticOptimizerParams = QuadraticOptimizerParams(),
                  ini_weights::Vector{Float64} = rand(nvars(problem)) # unused
)
    if nvars(problem) == 0
        return CoverProblemResult(Vector{Float64}(), Vector{Float64}(), 0.0, 0.0)
    end

    # Perform the optimization
    #try
    # using JuMP
    m = opt_model(problem)
    setsolver(m, quadratic_solver(opt_params))
    solve(m)
    w = copy(getvalue(getindex(m, :w)))
    # remove small non-zero probabilities due to optimization method errors
    for i in eachindex(w)
        @inbounds if w[i] < problem.params.min_weight
            w[i] = 0.0
        end
    end
    s = getobjectivevalue(m)
    return CoverProblemResult(w, problem.var_scores .* w, s, s)
    #catch x
    #    warn("Exception in optimize(CoverProblem): $x")
    #    return nothing
    #end
end

varXvar_mul!(vvXw::AbstractVector{Float64},
             problem::QuadraticCoverProblem,
             w::AbstractVector{Float64}) =
    A_mul_B!(vvXw, problem.varXvar_scores, w)
