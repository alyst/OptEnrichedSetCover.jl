const MultiobjCoverProblemFitnessScheme =
    ParetoFitnessScheme{2,Float64,true,MultiobjCoverProblemScoreAggregator{MultiobjProblemSoftFold2d}}

"""
Parameters for the Borg-based optimization of [`MultiobjCoverProblem`](@ref).

See [`optimize`](@ref).
"""
struct MultiobjOptimizerParams <: AbstractOptimizerParams{MultiobjCoverProblem}
    pop_size::Int
    weight_digits::Union{Int, Nothing}
    max_steps::Int
    max_steps_without_progress::Int
    fitness_tolerance::Float64
    min_delta_fitness_tolerance::Float64
    fold_ratio_threshold::Float64
    #workers::Vector{Int}
    nworkers::Int
    trace_interval::Float64

    borg_params::BlackBoxOptim.ParamsDict

    function MultiobjOptimizerParams(;
        # default Borg/BBO Opt.Controller parameter overrides
        PopulationSize::Integer = 100,
        WeightDigits::Union{Integer, Nothing} = 2, # precision digits for set weights
        MaxSteps::Integer = 10_000_000,
        MaxStepsWithoutProgress::Integer = 50_000,
        FitnessTolerance::Real = 0.1,
        MinDeltaFitnessTolerance::Real = 1E-5,
        FoldRatioThreshold::Real = 1.0,
        NWorkers::Integer = max(1, Threads.nthreads() - 1), #Workers::AbstractVector{Int} = Vector{Int}(),
        TraceInterval::Real = 5.0,
        kwargs...
    )
        if #=isempty(Workers) &&=# NWorkers > max(1, Threads.nthreads()-1)
            @warn("Requested NWorkers=$NWorkers, while only $(Threads.nthreads()) worker thread(s) available, reducing")
            NWorkers = Threads.nthreads()-1
        end
        #if isempty(Workers) && NWorkers > 1
        #    Workers = workers()[1:NWorkers]
        #end
        new(PopulationSize, WeightDigits,
            MaxSteps, MaxStepsWithoutProgress,
            FitnessTolerance, MinDeltaFitnessTolerance, FoldRatioThreshold,
            NWorkers, TraceInterval,
            BlackBoxOptim.kwargs2dict(kwargs...))
    end
end

function borg_params(opt_params::MultiobjOptimizerParams,
                     problem::MultiobjCoverProblem)
    if haskey(opt_params.borg_params, :ϵ)
        eps = opt_params.borg_params[:ϵ]
    else
        #eps_scale = get(opt_params.borg_params, :ϵ_scale, [0.2, 0.1, 0.1])
        #eps = eps_scale .* [score_scales(problem)...]
        eps = 0.01
    end
    BlackBoxOptim.chain(BlackBoxOptim.BorgMOEA_DefaultParameters,
        ParamsDict(:PopulationSize=>opt_params.pop_size,
                   :NThreads=>opt_params.nworkers > 1 ? opt_params.nworkers : 0,
                   :ϵ=>eps),
        opt_params.borg_params)
end

bbo_ctrl_params(opt_params::MultiobjOptimizerParams) =
    BlackBoxOptim.chain(BlackBoxOptim.DefaultParameters,
        ParamsDict(:PopulationSize=>opt_params.pop_size,
                   :MaxSteps=>opt_params.max_steps,
                   :MaxStepsWithoutProgress=>opt_params.max_steps_without_progress,
                   :FitnessTolerance=>opt_params.fitness_tolerance,
                   :MinDeltaFitnessTolerance=>opt_params.min_delta_fitness_tolerance,
                   :TraceInterval=>opt_params.trace_interval))

genop_params(opt_params::MultiobjOptimizerParams) =
    BlackBoxOptim.ParamsDict()

"""
Wraps [`MultiobjCoverProblem`](@ref) as [`BlackBoxOptim.OptimizationProblem`]
"""
struct MultiobjCoverProblemBBOWrapper{P <: MultiobjCoverProblem, SS <: RectSearchSpace} <:
        BlackBoxOptim.OptimizationProblem{MultiobjCoverProblemFitnessScheme}
    orig::P
    fitness_scheme::MultiobjCoverProblemFitnessScheme

    search_space::SS

    function MultiobjCoverProblemBBOWrapper(
        orig::MultiobjCoverProblem,
        score_folding::DetailedScoreFolding,
        dimdigits::Union{Integer,Nothing}=2
    )
        fitness_scheme = MultiobjCoverProblemFitnessScheme(
            aggregator=MultiobjCoverProblemScoreAggregator(orig.params, score_folding))
        ss = RectSearchSpace(nvars(orig), (0.0, 1.0), dimdigits=dimdigits)
        new{typeof(orig), typeof(ss)}(orig, fitness_scheme, ss)
    end
end

Base.copy(problem::MultiobjCoverProblemBBOWrapper) =
    MultiobjCoverProblemBBOWrapper(problem.orig, dimdigits(search_space(problem), 1))
#=
BlackBoxOptim.show_fitness(io::IO, score::DetailedScore,
                           problem::MultiobjCoverProblemBBOWrapper{MultiobjProblemNoFolding}) =
    @printf(io, "(sets=%.3f set×set=%.3f uncovered=%.3f covered=%.3f)\n",
            score[1], score[2], score[3], score[4])

BlackBoxOptim.show_fitness(io::IO, score::NTuple{2,Float64},
                           problem::MultiobjCoverProblemBBOWrapper{MultiobjProblemDrop3}) =
    @printf(io, "(sets=%.3f set×set=%.3f)\n", score[1], score[2]^2)

BlackBoxOptim.show_fitness(io::IO, score::NTuple{3,Float64},
                           problem::MultiobjCoverProblemBBOWrapper{MultiobjProblemFold3d}) =
    @printf(io, "(sets+k(sets²+uncov)=%.3f sets²+k(sets+uncov)=%.3f uncov+k(sets+sets²)=%.3f)\n",
          score[1], score[2], score[3])
=#

BlackBoxOptim.show_fitness(io::IO, score::FoldedScore,
                           problem::MultiobjCoverProblemBBOWrapper) =
    @printf(io, "(sets+k·miscover=%.3f sets²=%.3f)\n", score[1], score[2])

BlackBoxOptim.show_fitness(io::IO, score::IndexedTupleFitness, problem::MultiobjCoverProblemBBOWrapper) =
    BlackBoxOptim.show_fitness(io, score.orig, problem)

fold_score(score::DetailedScore, p::MultiobjCoverProblemBBOWrapper) =
    p.fitness_scheme.aggregator.score_folding(score)

BlackBoxOptim.fitness(x::BlackBoxOptim.AbstractIndividual,
                      p::MultiobjCoverProblemBBOWrapper) =
    fold_score(score(x, p.orig), p)

generate_recombinators(problem::MultiobjCoverProblemBBOWrapper, params) =
    CrossoverOperator[BlackBoxOptim.DiffEvoRandBin1(chain(BlackBoxOptim.DE_DefaultOptions, params)),
                    SimplexCrossover{3}(1.05),
                    SimplexCrossover{2}(1.1),
                    SimulatedBinaryCrossover(0.05, 16.0),
                    SimulatedBinaryCrossover(0.05, 3.0),
                    SimulatedBinaryCrossover(0.1, 5.0),
                    SimulatedBinaryCrossover(0.2, 16.0),
                    UnimodalNormalDistributionCrossover{2}(chain(BlackBoxOptim.UNDX_DefaultOptions, params)),
                    UnimodalNormalDistributionCrossover{3}(chain(BlackBoxOptim.UNDX_DefaultOptions, params)),
                    ParentCentricCrossover{2}(chain(BlackBoxOptim.PCX_DefaultOptions, params)),
                    ParentCentricCrossover{3}(chain(BlackBoxOptim.PCX_DefaultOptions, params)),
]

"""
Switch from one activated var to another one that is connected to it
by the set overlap.
"""
struct VarSwitchMutation <: BlackBoxOptim.MutationOperator
    varXvar_probs::Vector{Tuple{Vector{Int}, Weights{Float64}}} # var -> var switching probabilities
end

function VarSwitchMutation(problem::MultiobjCoverProblem;
                           score_shape::Real = 0.1, vXv_cutoff::Real = -5.0)
    score_shape > 0 || throw(ArgumentError("score_shape must be positive"))
    varXvar_probs = Vector{Tuple{Vector{Int}, Weights{Float64}}}()
    for j in axes(problem.varXvar_scores, 2)
        vixs = Vector{Int}()
        w = Vector{Float64}()
        for i in axes(problem.varXvar_scores, 1)
            i == j && continue
            @inbounds vXv = problem.varXvar_scores[i, j]
            @assert vXv <= 0
            if vXv < vXv_cutoff
                push!(vixs, i)
                push!(w, 1 - exp(vXv * score_shape))
            end
        end
        push!(varXvar_probs, (vixs, weights(w)))
    end
    return VarSwitchMutation(varXvar_probs)
end

function BlackBoxOptim.apply!(m::VarSwitchMutation, v::AbstractVector{Float64}, index::Int)
    nnz = length(v) - sum(iszero, v)
    (nnz == 0) && return v
    # select a random vertex ix with non-zero weight to switch
    nzix = sample(1:nnz)
    ix = findfirst(!iszero, v)
    for _ in 1:(nzix-1)
        ix = findnext(!iszero, v, ix)
    end
    # select the neighbouring vertex to switch to
    jixs, jws = m.varXvar_probs[ix]
    isempty(jixs) && return v # nowhere to switch to
    jx = sample(jixs, jws)
    # do the switch
    v[jx] = max(v[jx], v[ix])
    v[ix] = 0.0
    return v
end

generate_modifier(problem::MultiobjCoverProblemBBOWrapper, params) =
    FixedGeneticOperatorsMixture(GeneticOperator[
                                 VarSwitchMutation(problem.orig),
                                 MutationClock(PolynomialMutation(search_space(problem), 30.0), 1/numdims(problem)),
                                 MutationClock(UniformMutation(search_space(problem)), 1/numdims(problem))],
                                 [0.5, 0.3, 0.2])

"""
The result of [`optimize`](@ref).
Contains the solutions on the Pareto front: weights of the annotation terms and corresponding cover scores.
"""
struct MultiobjCoverProblemResult
    var2set::Vector{Int}        # indices of sets in the original SetMosaic
    varscores::Vector{Float64}  # scores of the sets

    varweights::Matrix{Float64} # weights of the vars for each solution on Pareto front
    scores::Vector{DetailedScore} # detailed scores for each solution
    folded_scores::Vector{FoldedScore}
    agg_scores::Vector{Float64}
    best_ix::Int                # index of the best solution (agg_score minimum)
end

function MultiobjCoverProblemResult(problem::MultiobjCoverProblem,
                                    optres::BlackBoxOptim.OptimizationResults)
    # initialize solutions
    wmtx = similar(best_candidate(optres), nvars(problem),
                   length(pareto_frontier(optres)))
    orig_folded_scores = [archived_fitness(sol).orig for sol in pareto_frontier(optres)]
    front_perm = sortperm(orig_folded_scores) # sort lexicographically
    @inbounds for i in axes(wmtx, 2)
        wmtx[:, i] = BlackBoxOptim.params(pareto_frontier(optres)[front_perm[i]])
    end
    # remove small non-zero probabilities due to optimization method errors
    minw = problem.params.min_weight
    maxw = 1.0 - problem.params.min_weight
    @inbounds for i in eachindex(wmtx)
        if wmtx[i] <= minw
            wmtx[i] = 0.0
        elseif wmtx[i] >= maxw
            wmtx[i] = 1.0
        end
    end
    # calculate detailed scores for the corrected weights, then fold and aggregate
    scores = dropdims(mapslices(w -> score(w, problem), wmtx; dims=1), dims=1)
    folded_scores = fitness_scheme(optres).aggregator.score_folding.(scores)
    agg_scores = aggscore.(scores, Ref(problem.params))

    return MultiobjCoverProblemResult(problem.var2set, problem.var_scores,
                                      wmtx, scores, folded_scores, agg_scores,
                                      findmin(agg_scores)[2])
end

nvars(res::MultiobjCoverProblemResult) = length(res.var2set)
nsolutions(res::MultiobjCoverProblemResult) = size(res.varweights, 2)

Base.isempty(res::MultiobjCoverProblemResult) = isempty(res.varweights)

varweights(res::MultiobjCoverProblemResult, sol_ix::Integer) =
    view(res.varweights, :, sol_ix)

best_index(res::MultiobjCoverProblemResult, ::Nothing=nothing) = res.best_ix

function best_index(res::MultiobjCoverProblemResult, params::CoverParams)
    # FIXME replace with findmin(aggscore, res) when Julia 1.0 support is dropped
    best_ix = 1
    min_score = aggscore(res.scores[best_ix], params)
    @inbounds for i in 2:nsolutions(res)
        sol_score = aggscore(res.scores[i], params)
        if sol_score < min_score
            best_ix = i
            min_score = sol_score
        end
    end
    return best_ix
end

best_varweights(result::MultiobjCoverProblemResult) =
    @inbounds(varweights(result, best_index(result)))

best_aggscore(result::MultiobjCoverProblemResult) =
    @inbounds(result.agg_scores[best_index(result)])

function set2var(result::MultiobjCoverProblemResult, setix::Int)
    varix = searchsortedlast(result.var2set, setix)
    if varix > 0 && result.var2set[varix] == setix
        return varix
    else
        return 0
    end
end

"""
    optimize(problem::MultiobjCoverProblem,
             [opt_params::MultiobjOptimizerParams]) -> MultiobjCoverProblemResult

Optimize [`MultiobjCoverProblem`](@ref) and return the result.
Uses Borf multi-objective optimization method from *BlackBoxOptim.jl* package.
"""
function optimize(problem::MultiobjCoverProblem,
                  opt_params::MultiobjOptimizerParams = MultiobjOptimizerParams())
    if nvars(problem) == 0
        return MultiobjCoverProblemResult(
            Vector{Int}(), Vector{Float64}(), Matrix{Float64}(undef, 0, 0),
            Vector{DetailedScore}(), Vector{FoldedScore}(), Vector{Float64}(), 0)
    end

    score_folding = MultiobjProblemSoftFold2d(problem.params, opt_params.fold_ratio_threshold)
    bbowrapper = MultiobjCoverProblemBBOWrapper(problem, score_folding, opt_params.weight_digits)
    popmatrix = rand_individuals(search_space(bbowrapper), opt_params.pop_size, method=:latin_hypercube)
    # two extreme solutions and one totally neutral
    size(popmatrix, 2) > 0 && (popmatrix[:, 1] .= 0.0)
    size(popmatrix, 2) > 1 && (popmatrix[:, 2] .= 1.0)
    size(popmatrix, 2) > 2 && (popmatrix[:, 3] .= 0.5)
    size(popmatrix, 2) > 3 && (popmatrix[:, 4] .= sortperm(problem.var_scores, rev=true)./length(problem.var_scores))
    # solutions that select topN most significant sets
    score_qtls = quantile(problem.var_scores, range(0.0, stop=1.0, length=10))
    for i in 1:10
        (size(popmatrix, 2) < i + 4) && break
        popmatrix[:, i+4] .= ifelse.(problem.var_scores .<= score_qtls[i], 1.0, 0.0)
    end
    N = fieldcount(fitness_type(fitness_scheme_type(score_folding))) # FIXME numobjectives(bbowrapper) gives wrong result, looks like Julia problem
    population = FitPopulation(popmatrix, nafitness(IndexedTupleFitness{N,Float64}), ntransient=1)

    go_params = genop_params(opt_params)
    bboptimizer = BlackBoxOptim.BorgMOEA(bbowrapper, population,
            generate_recombinators(bbowrapper, go_params),
            generate_modifier(bbowrapper, go_params),
            RandomBound(search_space(bbowrapper)),
            borg_params(opt_params, problem))

    bboctrl = BlackBoxOptim.OptController(bboptimizer, bbowrapper,
                 bbo_ctrl_params(opt_params))
    bbores = bboptimize(bboctrl)
    BlackBoxOptim.isinterrupted(bbores) && throw(InterruptException())
    return MultiobjCoverProblemResult(problem, bbores)
end
