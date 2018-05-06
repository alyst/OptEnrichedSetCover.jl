# defines how to fold fitness components
abstract type AbstractFitnessFolding{N} end

BlackBoxOptim.numobjectives(::Type{<:AbstractFitnessFolding{N}}) where N = N
BlackBoxOptim.numobjectives(::AbstractFitnessFolding{N}) where N = N

BlackBoxOptim.fitness_type(::Type{FF}) where FF <: AbstractFitnessFolding =
    fitness_type(fitness_scheme_type(FF))
BlackBoxOptim.fitness_type(::FF) where FF <: AbstractFitnessFolding =
    fitness_type(FF)
BlackBoxOptim.fitness_scheme_type(::FF) where FF <: AbstractFitnessFolding =
    fitness_scheme_type(FF)

abstract type MultiobjectiveProblemFitnessFolding{N} <: AbstractFitnessFolding{N} end

BlackBoxOptim.fitness_scheme_type(::Type{FF}) where FF <: MultiobjectiveProblemFitnessFolding{N} where N =
    ParetoFitnessScheme{N, Float64, true, MultiobjectiveCoverProblemScoreAggregator{FF}}

function MultiobjectiveProblemFitnessFolding(mosaic::MaskedSetMosaic, params::CoverParams,
                                             fitfolding::Union{Symbol, Void} = nothing)
    if fitfolding === nothing || fitfolding == :auto # chose folding automatically
        return MultiobjectiveProblemSoft12Convolute(params)
    elseif fitfolding == :none
        return MultiobjectiveProblemNoFolding()
    elseif fitfolding == :soft_convolute_sets_and_setXset
        return MultiobjectiveProblemSoft12Convolute(params)
    else
        throw(ArgumentError("Unknown fitfolding $fitfolding"))
    end
end

struct MultiobjectiveCoverProblemScoreAggregator{FF <: MultiobjectiveProblemFitnessFolding}
    k_sXs::Float64
    k_uncovered::Float64

    MultiobjectiveCoverProblemScoreAggregator{FF}(params::CoverParams) where FF =
        new{FF}(params.setXset_factor, NaN)
end

# no folding of fitness components
struct MultiobjectiveProblemNoFolding <: MultiobjectiveProblemFitnessFolding{3}
end
(::MultiobjectiveProblemNoFolding)(fitness::NTuple{2,Float64}) = fitness

(agg::MultiobjectiveCoverProblemScoreAggregator{MultiobjectiveProblemNoFolding})(score::NTuple{2, Float64}) =
    score[1] + agg.k_sXs * score[2]# + agg.k_mXm * score[3]^2

# sum var scores and setXset penalties (maskXmask penalties are separate)
struct MultiobjectiveProblemSoft12Convolute <: MultiobjectiveProblemFitnessFolding{2}
    setXset_factor::Float64

    MultiobjectiveProblemSoft12Convolute(params::CoverParams) =
        new(params.setXset_factor)
end

function (f::MultiobjectiveProblemSoft12Convolute)(fitness::NTuple{2,Float64})
    a = fitness[1]
    b = fitness[2]
    k = b/max(100.0, -a) - 3.0
    (k <= 0.0) && return (a, b)
    kk = 1.0 - 1.0/(1.0 + k)
    return (muladd(kk * f.setXset_factor, b, a), (1.0 - kk) * b)
end

(agg::MultiobjectiveCoverProblemScoreAggregator{MultiobjectiveProblemSoft12Convolute})(score::NTuple{2, Float64}) =
    score[1] + agg.k_sXs * score[2]

#=
# drop 3rd components (maskXmask)
struct MultiobjectiveProblemDrop3 <: MultiobjectiveProblemFitnessFolding{2}
end
(f::MultiobjectiveProblemDrop3)(fitness::NTuple{3,Float64}) =
    (fitness[1], sqrt(fitness[2]))
(agg::MultiobjectiveCoverProblemScoreAggregator{MultiobjectiveProblemDrop3})(score::NTuple{2, Float64}) =
    score[1] + agg.k_sXs * score[2]^2
=#

"""
Multi-objective optimal Enriched-Set Cover problem.
"""
struct MultiobjectiveCoverProblem{FF<:MultiobjectiveProblemFitnessFolding, F} <: AbstractCoverProblem{F}
    params::CoverParams
    fitfolding::FF

    var2set::Vector{Int}
    var_scores::Vector{Float64}
    varXvar_scores::Matrix{Float64}

    function MultiobjectiveCoverProblem(params::CoverParams,
                        fitfolding::FF,
                        var2set::AbstractVector{Int},
                        var_scores::AbstractVector{Float64},
                        varXvar_scores::AbstractMatrix{Float64}
    ) where {FF<:MultiobjectiveProblemFitnessFolding}
        length(var2set) == length(var_scores) ==
        size(varXvar_scores, 1) == size(varXvar_scores, 2) ||
            throw(ArgumentError("var2set, var_scores and varXvar_scores counts do not match"))
        F = fitness_type(FF)
        new{FF, F}(params, fitfolding, var2set, var_scores, varXvar_scores)
    end
end

nmasks(problem::MultiobjectiveCoverProblem) = length(problem.var_ranges)

BlackBoxOptim.fitness_scheme_type(::Type{<:MultiobjectiveCoverProblem{FF}}) where FF =
    BlackBoxOptim.fitness_scheme_type(FF)
BlackBoxOptim.fitness_type(::Type{<:MultiobjectiveCoverProblem{FF}}) where FF =
    BlackBoxOptim.fitness_type(FF)
BlackBoxOptim.numobjectives(::Type{<:MultiobjectiveCoverProblem{FF}}) where FF =
    numobjectives(FF)

function score_scales(problem::MultiobjectiveCoverProblem)
    w_min, w_max = extrema(problem.var_scores)
    sXs_min, sXs_max = length(problem.varXvar_scores) > 1 ? (Inf, -Inf) : (0.0, 0.0)
    return (max(w_max - w_min, 0.1),
            max(sXs_max - sXs_min, 0.1))
end

struct MultiobjectiveOptimizerParams <: AbstractOptimizerParams{MultiobjectiveCoverProblem}
    pop_size::Int
    max_steps::Int
    max_steps_without_progress::Int
    fitness_tolerance::Float64
    min_delta_fitness_tolerance::Float64
    trace_interval::Float64
    workers::Vector{Int}
    borg_params::BlackBoxOptim.ParamsDict

    function MultiobjectiveOptimizerParams(;
        # default Borg/BBO Opt.Controller parameter overrides
        NWorkers::Integer = 1, Workers::AbstractVector{Int} = Vector{Int}(),
        PopulationSize::Integer = 100,
        MaxSteps::Integer = 10_000_000,
        MaxStepsWithoutProgress::Integer = 10_000,
        FitnessTolerance::Real = 0.1,
        MinDeltaFitnessTolerance::Real = 1E-5,
        TraceInterval::Real = 5.0,
        kwargs...
    )
        if isempty(Workers) && NWorkers > 1
            Workers = workers()[1:NWorkers]
        end
        new(PopulationSize, MaxSteps, MaxStepsWithoutProgress,
            FitnessTolerance, MinDeltaFitnessTolerance, TraceInterval,
            Workers,
            BlackBoxOptim.kwargs2dict(kwargs))
    end
end

function borg_params(opt_params::MultiobjectiveOptimizerParams,
                     problem::MultiobjectiveCoverProblem)
    if haskey(opt_params.borg_params, :ϵ)
        eps = opt_params.borg_params[:ϵ]
    else
        #eps_scale = get(opt_params.borg_params, :ϵ_scale, [0.2, 0.1, 0.1])
        #eps = eps_scale .* [score_scales(problem)...]
        eps = [0.1, 0.1]
    end
    BlackBoxOptim.chain(BlackBoxOptim.BorgMOEA_DefaultParameters,
        ParamsDict(:PopulationSize=>opt_params.pop_size,
                   :Workers=>opt_params.workers,
                   :ϵ=>eps),
        opt_params.borg_params)
end

bbo_ctrl_params(opt_params::MultiobjectiveOptimizerParams) =
    BlackBoxOptim.chain(BlackBoxOptim.DefaultParameters,
        ParamsDict(:PopulationSize=>opt_params.pop_size,
                   :MaxSteps=>opt_params.max_steps,
                   :MaxStepsWithoutProgress=>opt_params.max_steps_without_progress,
                   :FitnessTolerance=>opt_params.fitness_tolerance,
                   :MinDeltaFitnessTolerance=>opt_params.min_delta_fitness_tolerance,
                   :TraceInterval=>opt_params.trace_interval))

genop_params(opt_params::MultiobjectiveOptimizerParams) =
    BlackBoxOptim.ParamsDict()

function MultiobjectiveCoverProblem(mosaic::MaskedSetMosaic, params::CoverParams = CoverParams();
                                    fitfolding::Union{Symbol, Void} = nothing)
    v2set = var2set(mosaic)
    v_scores = var_scores(mosaic, v2set, params)
    vXv_scores = varXvar_scores(mosaic, v2set, params, false)
    MultiobjectiveCoverProblem(params, MultiobjectiveProblemFitnessFolding(mosaic, params, fitfolding),
                               v2set, v_scores, vXv_scores)
end

# \sum_{i, j} A_{i, j} min(w_i, w_j)
function minplus_quad(A::AbstractMatrix{T},
                      w::AbstractVector{T}) where T
    length(w) == size(A, 1) == size(A, 2) ||
        throw(DimensionMismatch("A $(size(A)) and w ($(length(w))) size mismatch"))
    isempty(w) && return zero(T)
    res = zero(T)
    @inbounds for i in eachindex(w)
        const wi = w[i]
        const ioffset = (i-1)*size(A, 1)
        resi = zero(T)
        if zero(T) < wi < one(T)
            @simd for j in eachindex(w)
                resi += A[ioffset+j]*min(wi, w[j])
            end
        elseif wi == one(T) # wi >= maxw
            @simd for j in eachindex(w)
                resi += A[ioffset+j] * w[j] # FIXME dotprod() faster?
            end
        end
        res += resi
    end
    return res
end

# calculate r_i = fold(r_i, \sum_{j = 1} A_{i, j} min(u_i, v_j)), i = 1..N
function minplus_bilinear!(foldl::Function,
                           res::AbstractVector{T},
                           A::AbstractMatrix{T},
                           u::AbstractVector{T},
                           v::AbstractVector{T}) where T
    length(res) == length(u) == size(A, 2) ||
        throw(DimensionMismatch("res ($(length(res))), u ($(length(u))) or A $(size(A)) size mismatch"))
    length(v) == size(A, 1) ||
        throw(DimensionMismatch("v ($(length(v))) and A $(size(A)) size mismatch"))
    n = length(res)
    @inbounds for i in 1:size(A, 2)
        const ui = u[i]
        const ioffset = (i-1)*size(A, 1)
        x = zero(T)
        if zero(T) < ui < one(T)
            @simd for j in eachindex(v)
                x += A[ioffset+j] * min(ui, v[j])
            end
        elseif ui == one(T)
            @simd for j in eachindex(v)
                x += A[ioffset+j] * v[j] # FIXME dotprod() faster?
            end
        end
        res[i] = foldl(res[i], x)
    end
    return res
end

take2(u, v) = v

function exclude_vars(problem::MultiobjectiveCoverProblem,
                      vars::AbstractVector{Int};
                      penalize_overlaps::Bool = true)
    varmask = fill(true, nvars(problem))
    varmask[vars] = false
    v_scores = problem.var_scores[varmask]
    if penalize_overlaps
        # penalize overlapping sets
        pweights = fill(0.0, nvars(problem))
        pweights[vars] = 1.0
        varscore_penalties = similar(pweights)
        minplus_bilinear!(take2, varscore_penalties,
                          problem.varXvar_scores,
                          pweights, pweights)
        v_scores .-= view(varscore_penalties, varmask) .* problem.params.setXset_factor
    end
    return MultiobjectiveCoverProblem(problem.params, problem.fitfolding,
                problem.var2set[varmask],
                v_scores, problem.varXvar_scores[varmask, varmask])
end

"""
Unfolded multiobjective score (fitness) of the OESC coverage.

* `w` probabilities of the sets being covered
"""
function rawscore(problem::MultiobjectiveCoverProblem, w::AbstractVector{Float64})
    @assert length(w) == nvars(problem)
    a = dot(problem.var_scores, w)
    if problem.params.setXset_factor == 0.0
        b = 0.0
    else
        # skip set interactions as it would be aggregated to zero
        setXset = fill!(similar(w), 0.0)
        minplus_bilinear!(min, setXset, problem.varXvar_scores, w, w)
        b = -sum(setXset)
    end
    return (a, b)
end

"""
Score (fitness) of the OESC coverage.

* `w` probabilities of the sets being covered
"""
score(problem::MultiobjectiveCoverProblem, w::AbstractVector{Float64}) =
    problem.fitfolding(rawscore(problem, w))

aggscore(problem::MultiobjectiveCoverProblem{FF}, w::AbstractVector{Float64}) where FF =
    MultiobjectiveCoverProblemScoreAggregator{FF}(problem.params)(score(problem, w))

# wraps MultiobjectiveCoverProblem as BlackBoxOptim OptimizationProblem
struct MultiobjectiveCoverProblemBBOWrapper{FF <: MultiobjectiveProblemFitnessFolding, FS <: FitnessScheme, F} <:
        BlackBoxOptim.OptimizationProblem{FS}
    orig::MultiobjectiveCoverProblem{FF, F}
    fitness_scheme::FS
    search_space::SearchSpace

    function MultiobjectiveCoverProblemBBOWrapper(
            orig::MultiobjectiveCoverProblem{FF, F}; digits::Integer=2
    ) where {FF, F}
        FS = fitness_scheme_type(FF)
        fitscheme = FS(aggregator=MultiobjectiveCoverProblemScoreAggregator{FF}(orig.params))
        new{FF, FS, F}(orig, fitscheme, symmetric_search_space(nvars(orig), (0.0, 1.0), digits=digits))
    end
end

Base.copy(problem::MultiobjectiveCoverProblemBBOWrapper) =
    MultiobjectiveCoverProblemBBOWrapper(problem.orig; digits=first(digits(search_space(problem))))

BlackBoxOptim.show_fitness(io::IO, score::NTuple{2,Float64},
                           problem::MultiobjectiveCoverProblemBBOWrapper{MultiobjectiveProblemNoFolding}) =
    @printf(io, "(sets=%.3f set×set=%.3f)\n", score[1], score[2])

#=
BlackBoxOptim.show_fitness(io::IO, score::NTuple{2,Float64},
                           problem::MultiobjectiveCoverProblemBBOWrapper{MultiobjectiveProblemDrop3}) =
    @printf(io, "(sets=%.3f set×set=%.3f)\n", score[1], score[2]^2)
=#
BlackBoxOptim.show_fitness(io::IO, score::NTuple{2,Float64},
                           problem::MultiobjectiveCoverProblemBBOWrapper{MultiobjectiveProblemSoft12Convolute}) =
  @printf(io, "(sets+k⋅set×set=%.3f (1.0-k)⋅set×set=%.3f)\n", score[1], score[2])

BlackBoxOptim.show_fitness(io::IO, score::IndexedTupleFitness, problem::MultiobjectiveCoverProblemBBOWrapper) =
    BlackBoxOptim.show_fitness(io, score.orig, problem)

BlackBoxOptim.fitness(x::BlackBoxOptim.AbstractIndividual,
                      p::MultiobjectiveCoverProblemBBOWrapper) = score(p.orig, x)

generate_recombinators(problem::OptimizationProblem, params) =
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

generate_modifier(problem::OptimizationProblem, params) =
    FixedGeneticOperatorsMixture(GeneticOperator[
                                 MutationClock(PolynomialMutation(search_space(problem), 30.0), 1/numdims(problem)),
                                 MutationClock(UniformMutation(search_space(problem)), 1/numdims(problem))],
                                 [0.75, 0.25])

function optimize(problem::MultiobjectiveCoverProblem,
                  opt_params::MultiobjectiveOptimizerParams = MultiobjectiveOptimizerParams())
    if nvars(problem) == 0
        w = Vector{Float64}()
        return CoverProblemResult(Vector{Int}(), w, Vector{Float64}(),
                                  score(problem, w), 0.0, nothing)
    end

    bbowrapper = MultiobjectiveCoverProblemBBOWrapper(problem)
    popmatrix = BlackBoxOptim.rand_individuals_lhs(search_space(bbowrapper), opt_params.pop_size)
    # two extreme solutions and one totally neutral
    size(popmatrix, 2) > 0 && (popmatrix[:, 1] = 0.0)
    size(popmatrix, 2) > 1 && (popmatrix[:, 2] = 1.0)
    size(popmatrix, 2) > 2 && (popmatrix[:, 3] = 0.5)
    size(popmatrix, 2) > 3 && (popmatrix[:, 4] = sortperm(problem.var_scores, rev=true)./length(problem.var_scores))
    # solutions that select topN most significant sets
    score_qtls = quantile(problem.var_scores, linspace(0.0, 1.0, 10))
    for i in 1:10
        (size(popmatrix, 2) < i + 4) && break
        popmatrix[:, i+4] = ifelse.(problem.var_scores .<= score_qtls[i], 1.0, 0.0)
    end
    const N = numobjectives(typeof(problem))
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
    isinterrupted(bbores) && throw(InterruptException())
    w = best_candidate(bbores)
    s = best_fitness(bbores)

    # remove small non-zero probabilities due to optimization method errors
    const minw = problem.params.min_weight
    const maxw = 1.0 - problem.params.min_weight
    @inbounds for i in eachindex(w)
        if w[i] <= minw
            w[i] = 0.0
        elseif w[i] >= maxw
            w[i] = 1.0
        end
    end

    fitness_frontier = [af.orig for af in archived_fitness.(pareto_frontier(bbores))]
    raw_fitness_frontier = [rawscore(problem, w) for w in BlackBoxOptim.params.(pareto_frontier(bbores))]
    frontier_perm = sortperm(raw_fitness_frontier)

    return CoverProblemResult(problem.var2set, w, problem.var_scores .* w, s,
                              BlackBoxOptim.aggregate(s, fitness_scheme(bbores)),
                              (fitness_frontier[frontier_perm],
                               raw_fitness_frontier[frontier_perm]))
end
