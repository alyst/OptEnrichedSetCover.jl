# defines how to fold fitness into N components
abstract type FitnessFolding{N} end

BlackBoxOptim.fitness_scheme_type(fitfold::FitnessFolding) =
    fitness_scheme_type(typeof(fitfold))

BlackBoxOptim.fitness_scheme_type(::Type{FF}) where FF <: FitnessFolding{N} where N =
    ParetoFitnessScheme{N, Float64, true, MultiobjCoverProblemScoreAggregator{FF}}

function MultiobjProblemFitnessFolding(params::CoverParams,
                                       kind::Union{Symbol, Nothing} = nothing)
    if kind === nothing || kind == :auto # chose folding automatically
        return MultiobjProblemSoftFold2d(params)
    elseif kind == :none
        return MultiobjProblemNoFolding()
    elseif kind == :fold2d
        return MultiobjProblemSoftFold2d(params)
    elseif kind == :fold3d
        return MultiobjProblemSoftFold3d(params)
    else
        throw(ArgumentError("Unknown fitness folding kind $kind"))
    end
end

struct MultiobjCoverProblemScoreAggregator{FF <: FitnessFolding}
    setXset_factor::Float64
    uncovered_factor::Float64
    covered_factor::Float64

    MultiobjCoverProblemScoreAggregator{FF}(params::CoverParams) where FF =
        new{FF}(params.setXset_factor,
                params.uncovered_factor,
                params.covered_factor)
end

# no folding of fitness components
struct MultiobjProblemNoFolding <: FitnessFolding{4}
end
(::MultiobjProblemNoFolding)(fitness::NTuple{4,Float64}) = fitness

(agg::MultiobjCoverProblemScoreAggregator{MultiobjProblemNoFolding})(score::NTuple{4, Float64}) =
    score[1] + agg.setXset_factor * score[2] +
    agg.uncovered_factor * score[3] + agg.covered_factor * score[4]

# sum var scores and setXset penalties (maskXmask penalties are separate)
struct MultiobjProblemSoftFold2d <: FitnessFolding{2}
    k::Float64
    k_sXs::Float64
    setXset_factor::Float64
    uncovered_factor::Float64
    covered_factor::Float64

    MultiobjProblemSoftFold2d(params::CoverParams, k::Float64=0.1) =
        new(params.setXset_factor > 0.0 ? k : 0.0,
            params.setXset_factor > 0.0 ? k/params.setXset_factor : 0.0,
            params.setXset_factor,
            params.setXset_factor > 0.0 ? params.uncovered_factor/params.setXset_factor : 1.0,
            params.setXset_factor > 0.0 ? params.covered_factor/params.setXset_factor : 1.0)
end

function (f::MultiobjProblemSoftFold2d)(fitness::NTuple{4,Float64})
    a = fitness[1]
    b = fitness[2] +
        f.uncovered_factor * fitness[3] +
        f.covered_factor * fitness[4]
    return (f.k * f.setXset_factor * b + (1-f.k) * a,
            (1-f.k) * b + f.k_sXs * a)
end

(agg::MultiobjCoverProblemScoreAggregator{MultiobjProblemSoftFold2d})(score::NTuple{2, Float64}) =
    score[1] + agg.setXset_factor * score[2]

# sum var scores and setXset penalties (maskXmask penalties are separate)
struct MultiobjProblemSoftFold3d <: FitnessFolding{3}
    k::Float64
    k_sXs::Float64
    k_u::Float64
    setXset_factor::Float64
    uncovered_factor::Float64
    covered_factor::Float64

    MultiobjProblemSoftFold3d(params::CoverParams, k::Float64=0.1) =
        new(params.setXset_factor > 0.0 ? k : 0.0,
            params.setXset_factor > 0.0 ? k/params.setXset_factor : 0.0,
            params.uncovered_factor > 0.0 ? k/params.uncovered_factor : 0.0,
            params.setXset_factor,
            params.uncovered_factor,
            params.uncovered_factor > 0.0 ? params.covered_factor/params.uncovered_factor : 1.0)
end

function (f::MultiobjProblemSoftFold3d)(fitness::NTuple{4,Float64})
    a = fitness[1]
    b = fitness[2]
    c = fitness[3] + f.covered_factor * fitness[4]
    return (0.5 * f.k * (f.setXset_factor * b + f.uncovered_factor * c) + (1-f.k) * a,
            0.5 * f.k_sXs * (a + f.uncovered_factor * c) + (1-f.k) * b,
            0.5 * f.k_u * (a + f.setXset_factor * b) + (1-f.k) * c)
end

(agg::MultiobjCoverProblemScoreAggregator{MultiobjProblemSoftFold3d})(score::NTuple{3, Float64}) =
    score[1] + agg.setXset_factor * score[2] + agg.uncovered_factor * score[3]

"""
Multi-objective optimal Enriched-Set Cover problem.
"""
struct MultiobjCoverProblem <: AbstractCoverProblem{NTuple{4, Float64}}
    params::CoverParams

    var2set::Vector{Int}

    var_scores::Vector{Float64}
    varXvar_scores::Matrix{Float64}

    tileXvar::SparseMaskMatrix      # tile-in-mask X set(var)
    nmasked_pertile::Vector{Int}    # number of masked elements for each tile
    nunmasked_pertile::Vector{Int}  # number of unmasked elements for each tile

    tilepool::Vector{Vector{Float64}}   # pool of tile-sized vectors FIXME move to evaluator
    varpool::Vector{Vector{Float64}}    # pool of var-sized vectors FIXME move to evaluator

    function MultiobjCoverProblem(params::CoverParams,
                        var2set::AbstractVector{Int},
                        var_scores::AbstractVector{Float64},
                        varXvar_scores::AbstractMatrix{Float64},
                        tileXvar::SparseMaskMatrix,
                        nmasked_pertile::AbstractVector{Int},
                        nunmasked_pertile::AbstractVector{Int}
    )
        length(var2set) == length(var_scores) ==
        size(varXvar_scores, 1) == size(varXvar_scores, 2) ==
        size(tileXvar, 2) ||
            throw(ArgumentError("var2set, var_scores and varXvar_scores counts do not match"))
        size(tileXvar, 1) == length(nmasked_pertile) == length(nunmasked_pertile) ||
            throw(ArgumentError("tileXvar and nmasked_pertile rows count does not match"))
        new(params, var2set,
            var_scores, varXvar_scores,
            tileXvar, nmasked_pertile, nunmasked_pertile,
            Vector{Vector{Float64}}(), Vector{Vector{Float64}}())
    end
end

ntiles(problem::MultiobjCoverProblem) = length(problem.nmasked_pertile)

function score_scales(problem::MultiobjCoverProblem)
    w_min, w_max = extrema(problem.var_scores)
    sXs_min, sXs_max = length(problem.varXvar_scores) > 1 ? (Inf, -Inf) : (0.0, 0.0)
    return (max(w_max - w_min, 0.1),
            max(sXs_max - sXs_min, 0.1))
end

struct MultiobjOptimizerParams <: AbstractOptimizerParams{MultiobjCoverProblem}
    fitness_folding::Union{Symbol, Nothing}
    pop_size::Int
    weight_digits::Union{Int, Nothing}
    max_steps::Int
    max_steps_without_progress::Int
    fitness_tolerance::Float64
    min_delta_fitness_tolerance::Float64
    trace_interval::Float64
    workers::Vector{Int}
    borg_params::BlackBoxOptim.ParamsDict

    function MultiobjOptimizerParams(;
        FitnessFolding::Union{Symbol, Nothing} = nothing,
        # default Borg/BBO Opt.Controller parameter overrides
        NWorkers::Integer = 1, Workers::AbstractVector{Int} = Vector{Int}(),
        PopulationSize::Integer = 100,
        WeightDigits::Union{Integer, Nothing} = 2, # precision digits for set weights
        MaxSteps::Integer = 10_000_000,
        MaxStepsWithoutProgress::Integer = 50_000,
        FitnessTolerance::Real = 0.1,
        MinDeltaFitnessTolerance::Real = 1E-5,
        TraceInterval::Real = 5.0,
        kwargs...
    )
        if isempty(Workers) && NWorkers > nworkers()
            @warn("Requested NWorkers=$NWorkers, while only $(nworkers()) available, reducing")
            NWorkers = nworkers()
        end
        if isempty(Workers) && NWorkers > 1
            Workers = workers()[1:NWorkers]
        end
        new(FitnessFolding, PopulationSize, WeightDigits, MaxSteps, MaxStepsWithoutProgress,
            FitnessTolerance, MinDeltaFitnessTolerance, TraceInterval,
            Workers,
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
                   :Workers=>opt_params.workers,
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

function MultiobjCoverProblem(mosaic::MaskedSetMosaic, params::CoverParams = CoverParams())
    v2set = var2set(mosaic)
    v_scores = var_scores(mosaic, v2set, params)
    vXv_scores = varXvar_scores(mosaic, v2set, params, false)
    tXv, nmasked_pertile, nunmasked_pertile = tilemaskXvar(mosaic)
    MultiobjCoverProblem(params,
                         v2set, v_scores, vXv_scores,
                         tXv, nmasked_pertile, nunmasked_pertile)
end

# \sum_{i, j} A_{i, j} min(w_i, w_j)
function minplus_quad(A::AbstractMatrix{T},
                      w::AbstractVector{T}) where T
    length(w) == size(A, 1) == size(A, 2) ||
        throw(DimensionMismatch("A $(size(A)) and w ($(length(w))) size mismatch"))
    isempty(w) && return zero(T)
    res = zero(T)
    @inbounds for i in eachindex(w)
        wi = w[i]
        ioffset = (i-1)*size(A, 1)
        resi = zero(T)
        if zero(T) < wi < one(T)
            for j in eachindex(w)
                resi += A[ioffset+j]*min(wi, w[j])
            end
        elseif wi == one(T) # wi >= maxw
            for j in eachindex(w)
                resi += A[ioffset+j] * w[j] # FIXME dotprod() faster?
            end
        end
        res += resi
    end
    return res
end

# calculate res_i = fold(res_i, \sum_{j = 1} A_{i, j} min(u_i, v_j)), i = 1..N
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
        ui = u[i]
        ioffset = (i-1)*size(A, 1)
        x = zero(T)
        if zero(T) < ui < one(T) # if ui == 0, nothing is done
            for j in eachindex(v)
                x += A[ioffset+j] * min(ui, v[j])
            end
        elseif ui == one(T)
            for j in eachindex(v)
                x += A[ioffset+j] * v[j] # FIXME dotprod() faster?
            end
        end
        res[i] = foldl(res[i], x)
    end
    return res
end

take2(u, v) = v

function exclude_vars(problem::MultiobjCoverProblem,
                      vars::AbstractVector{Int};
                      penalize_overlaps::Bool = true)
    varmask = fill(true, nvars(problem))
    varmask[vars] .= false
    v_scores = problem.var_scores[varmask]
    if penalize_overlaps
        # penalize overlapping sets
        pweights = fill(0.0, nvars(problem))
        pweights[vars] .= 1.0
        varscore_penalties = similar(pweights)
        minplus_bilinear!(take2, varscore_penalties,
                          problem.varXvar_scores,
                          pweights, pweights)
        v_scores .-= view(varscore_penalties, varmask) .* problem.params.setXset_factor
    end
    return MultiobjCoverProblem(problem.params,
                problem.var2set[varmask],
                v_scores, problem.varXvar_scores[varmask, varmask],
                problem.tileXvar[:, varmask], # FIXME can also remove unused tiles
                problem.nmasked_pertile, problem.nunmasked_pertile)
end

# the tuple of:
# - total number of uncovered masked elements (in all masks overlapping with the cover)
# -`` total number of covered unmasked elements (in all masks overlapping with the cover)
function miscover_score(problem::MultiobjCoverProblem, w::AbstractVector{Float64})
    __check_vars(w, problem)
    isempty(problem.nmasked_pertile) && return (0.0, 0.0)
    # calculate the tile coverage weights
    wtile = isempty(problem.tilepool) ?
            Vector{Float64}(undef, ntiles(problem)) :
            pop!(problem.tilepool)
    fill!(wtile, 0.0)
    @inbounds for (varix, wvar) in enumerate(w)
        if wvar == 1.0
            wtile[view(problem.tileXvar, :, varix)] .= 1.0
        elseif wvar > 0.0
            for tileix in view(problem.tileXvar, :, varix)
                wtile[tileix] = max(wtile[tileix], wvar)
            end
        end
    end
    res = sum(problem.nmasked_pertile) - dot(problem.nmasked_pertile, wtile),
          dot(problem.nunmasked_pertile, wtile)
    push!(problem.tilepool, wtile) # release to the pool
    return res
end

"""
Unfolded multiobjective score (fitness) of the OESC coverage.

* `w` probabilities of the sets being covered
"""
function score(problem::MultiobjCoverProblem, w::AbstractVector{Float64})
    __check_vars(w, problem)
    a = dot(problem.var_scores, w)
    if problem.params.setXset_factor == 0.0
        b = 0.0
    else
        # skip set interactions as it would be aggregated to zero
        setXset = isempty(problem.varpool) ? similar(w) : pop!(problem.varpool)
        fill!(setXset, 0.0)
        minplus_bilinear!(min, setXset, problem.varXvar_scores, w, w)
        b = -sum(setXset)
        push!(problem.varpool, setXset)
    end
    if problem.params.covered_factor == 0.0 &&
       problem.params.uncovered_factor == 0.0
        # skip miscover score calculation if it's excluded from aggregation
        c = 0.0
        d = 0.0
    else
        c, d = miscover_score(problem, w)
    end
    return (a, b, c, d)
end

function aggscore(problem::MultiobjCoverProblem, w::AbstractVector{Float64})
    s = score(problem, w)
    s[1] +
    problem.params.setXset_factor * s[2] +
    problem.params.uncovered_factor * s[3] +
    problem.params.covered_factor * s[4]
end

# wraps MultiobjCoverProblem as BlackBoxOptim OptimizationProblem
struct MultiobjCoverProblemBBOWrapper{FF <: FitnessFolding, FS <: FitnessScheme, SS <: SearchSpace} <:
        BlackBoxOptim.OptimizationProblem{FS}
    orig::MultiobjCoverProblem
    fitness_folding::FF
    fitness_scheme::FS
    search_space::SS

    function MultiobjCoverProblemBBOWrapper(
        orig::MultiobjCoverProblem,
        fitness_folding::FF,
        dimdigits::Union{Integer,Nothing}=2
    ) where FF
        FS = fitness_scheme_type(FF)
        fitness_scheme = FS(aggregator=MultiobjCoverProblemScoreAggregator{FF}(orig.params))
        ss = RectSearchSpace(nvars(orig), (0.0, 1.0), dimdigits=dimdigits)
        new{FF, FS, typeof(ss)}(orig, fitness_folding, fitness_scheme, ss)
    end
end

Base.copy(problem::MultiobjCoverProblemBBOWrapper) =
    MultiobjCoverProblemBBOWrapper(problem.orig, problem.fitness_folding,
                                   dimdigits(search_space(problem), 1))

BlackBoxOptim.show_fitness(io::IO, score::NTuple{4,Float64},
                           problem::MultiobjCoverProblemBBOWrapper{MultiobjProblemNoFolding}) =
    @printf(io, "(sets=%.3f set×set=%.3f uncovered=%.3f covered=%.3f)\n",
            score[1], score[2], score[3], score[4])

#=
BlackBoxOptim.show_fitness(io::IO, score::NTuple{2,Float64},
                           problem::MultiobjCoverProblemBBOWrapper{MultiobjProblemDrop3}) =
    @printf(io, "(sets=%.3f set×set=%.3f)\n", score[1], score[2]^2)
=#
BlackBoxOptim.show_fitness(io::IO, score::NTuple{2,Float64},
                           problem::MultiobjCoverProblemBBOWrapper{MultiobjProblemSoftFold2d}) =
    @printf(io, "(sets+k·set²=%.3f set²+k·sets=%.3f)\n", score[1], score[2])

BlackBoxOptim.show_fitness(io::IO, score::NTuple{3,Float64},
                           problem::MultiobjCoverProblemBBOWrapper{MultiobjProblemSoftFold3d}) =
    @printf(io, "(sets+k(set²+uncov)=%.3f set²+k(sets+uncov)=%.3f uncov+k(sets+set²)=%.3f)\n",
          score[1], score[2], score[3])

BlackBoxOptim.show_fitness(io::IO, score::IndexedTupleFitness, problem::MultiobjCoverProblemBBOWrapper) =
    BlackBoxOptim.show_fitness(io, score.orig, problem)

BlackBoxOptim.fitness(x::BlackBoxOptim.AbstractIndividual,
                      p::MultiobjCoverProblemBBOWrapper) = p.fitness_folding(score(p.orig, x))

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

function optimize(problem::MultiobjCoverProblem,
                  opt_params::MultiobjOptimizerParams = MultiobjOptimizerParams())
    if nvars(problem) == 0
        w = Vector{Float64}()
        return CoverProblemResult(Vector{Int}(), w, Vector{Float64}(),
                                  score(problem, w), 0.0, nothing)
    end

    fitfolding = MultiobjProblemFitnessFolding(problem.params, opt_params.fitness_folding)
    bbowrapper = MultiobjCoverProblemBBOWrapper(problem, fitfolding,
                                                opt_params.weight_digits)
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
    N = fieldcount(fitness_type(fitness_scheme_type(fitfolding))) # FIXME numobjectives(bbowrapper) gives wrong result, looks like Julia problem
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
    w = best_candidate(bbores)

    # remove small non-zero probabilities due to optimization method errors
    minw = problem.params.min_weight
    maxw = 1.0 - problem.params.min_weight
    @inbounds for i in eachindex(w)
        if w[i] <= minw
            w[i] = 0.0
        elseif w[i] >= maxw
            w[i] = 1.0
        end
    end

    fitness_frontier = [archived_fitness(af).orig for af in pareto_frontier(bbores)]
    raw_fitness_frontier = score.(Ref(problem), BlackBoxOptim.params.(pareto_frontier(bbores)))
    frontier_perm = sortperm(raw_fitness_frontier)

    return CoverProblemResult(problem.var2set, w, problem.var_scores .* w,
                              score(problem, w), aggscore(problem, w),
                              (fitness_frontier[frontier_perm],
                               raw_fitness_frontier[frontier_perm]))
end
