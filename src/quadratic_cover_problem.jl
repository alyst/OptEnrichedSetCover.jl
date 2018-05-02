"""
Quadratic Programming Optimal Enriched-Set Cover problem.
"""
struct QuadraticCoverProblem <: AbstractCoverProblem{Float64}
    params::CoverParams

    var2set::Vector{Int}
    var_scores::Vector{Float64}
    varXvar_scores::Matrix{Float64}

    function QuadraticCoverProblem(params::CoverParams,
            var2set::Vector{Int},
            var_scores::Vector{Float64},
            varXvar_scores::Matrix{Float64}
    )
        length(var2set) == length(var_scores) ==
        size(varXvar_scores, 1) == size(varXvar_scores, 2) ||
            throw(ArgumentError("var2set, var_scores and varXvar_scores sizes do not match"))
        new(params, var2set, var_scores, varXvar_scores)
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
    v2set = var2set(mosaic)
    QuadraticCoverProblem(params, v2set,
                          var_scores(mosaic, v2set, params),
                          varXvar_scores(mosaic, v2set, params, true))
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
        return CoverProblemResult(problem.var2set, Vector{Float64}(), Vector{Float64}(), 0.0, 0.0)
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
    return CoverProblemResult(problem.var2set, w, problem.var_scores .* w, s, s)
    #catch x
    #    warn("Exception in optimize(CoverProblem): $x")
    #    return nothing
    #end
end

function exclude_vars(problem::QuadraticCoverProblem,
                      vars::AbstractVector{Int};
                      penalize_overlaps::Bool = true)
    varmask = fill(true, nvars(problem))
    varmask[vars] = false
    v_scores = problem.var_scores[varmask]
    if penalize_overlaps
        # penalize overlapping sets
        penalty_weights = fill(0.0, nvars(problem))
        penalty_weights[vars] = 1.0
        v_scores .-= (problem.varXvar_scores * penalty_weights)[varmask]
    end
    return QuadraticCoverProblem(problem.params,
                problem.var2set[varmask],
                v_scores, problem.varXvar_scores[varmask, varmask])
end
