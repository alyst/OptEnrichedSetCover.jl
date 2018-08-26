# Optimal Enriched-Set Cover
module OptEnrichedSetCover

using Requires, DataFrames, Distributions, BlackBoxOptim
using LinearAlgebra, Distributed
using Printf: @printf

export SetMosaic, CoverParams,
    AbstractCoverProblem, QuadraticCoverProblem, MultiobjectiveCoverProblem,
    AbstractOptimizerParams, QuadraticOptimizerParams, MultiobjectiveOptimizerParams,
    CoverEnumerationParams, CoverCollection, CoverProblemResult,
    nelements, ntiles, nsets, nmaskedsets, nmasks, nvars,
    nmasked, nunmasked, maskedset,
    tile, tiles, set, setsize,
    unmask, mask,
    aggscore, score, # any conflicts
    optimize # conflicts with Optim.jl

global __quadratic_problem_supported__ = false
global __ipopt_solver_available__ = false
global with_default_quadratic_solver = () -> error("Default quadratic solver not defined")

function ipopt_quadratic_solver()
    __quadratic_problem_supported__ ||
        error("Quadratic problem is not supported (JuMP or MathProgBase packages missing)")
    __ipopt_solver_available__ ||
        error("Ipopt package is missing")
    with_optimizer(IpoptSolver, print_level=0)
end

#using Mosek
#with_default_quadratic_solver() = with_optimizer(MosekSolver)
#using Gurobi
#with_default_quadratic_solver() = with_optimizer(GurobiSolver, OutputFlag=0)
#using CPLEX
#with_default_quadratic_solver() = with_optimizer(CplexSolver)

import BlackBoxOptim: OptimizationProblem

include("sparse_mask_matrix.jl")
include("set_score.jl")
include("mosaic.jl")
include("masked_mosaic.jl")
include("cover_problem.jl")
include("quadratic_cover_problem.jl")
include("multiobj_cover_problem.jl")
include("cover_enumerator.jl")
include("parallel.jl")

function __init__()
    # initialize optional dependencies
    has_JuMP = false
    @require JuMP="4076af6c-e467-56ae-b986-b466b2749572" begin
        using JuMP
        has_JuMP = true
    end
    @require Ipopt="b6b21f68-93f8-5de0-b562-5493be1d77c9" begin
        using Ipopt
        global __ipopt_solver_available__ = true
    end
    has_MathProgBase = false
    @require MathProgBase="fdba3010-5040-5b88-9595-932c9decdf73" begin
        using MathProgBase
        has_MathProgBase = true
    end

    global __quadratic_problem_supported__ = has_JuMP && has_MathProgBase && __ipopt_solver_available__
    # select the default solver for quadratic problem
    if __quadratic_problem_supported__ && __ipopt_solver_available__
        global with_default_quadratic_solver = ipopt_quadratic_solver
    end
end

end # module
