# defines how to fold fitness into N components
abstract type FitnessFolding{N} end

# raw cover score:
# 1 => covered set scores,
# 2 => covered setXset penalties,
# 3 => uncovered masked elements penalties
# 4 => covered unmasked elements penalties
const RawScore = NTuple{4, Float64}

# soft Fold2D score
const FoldedScore = NTuple{2, Float64}

BlackBoxOptim.fitness_scheme_type(fitfold::FitnessFolding) =
    fitness_scheme_type(typeof(fitfold))

BlackBoxOptim.fitness_scheme_type(::Type{FF}) where FF <: FitnessFolding{N} where N =
    ParetoFitnessScheme{N, Float64, true, MultiobjCoverProblemScoreAggregator{FF}}

struct MultiobjCoverProblemScoreAggregator{FF <: FitnessFolding}
    setXset_factor::Float64
    uncovered_factor::Float64
    covered_factor::Float64

    MultiobjCoverProblemScoreAggregator{FF}(params::CoverParams) where FF =
        new{FF}(params.setXset_factor,
                params.uncovered_factor,
                params.covered_factor)
end

#=
# no folding of fitness components
struct MultiobjProblemNoFolding <: FitnessFolding{4}
end
(::MultiobjProblemNoFolding)(fitness::RawScore) = fitness

(agg::MultiobjCoverProblemScoreAggregator{MultiobjProblemNoFolding})(score::RawScore) =
    score[1] + agg.setXset_factor * score[2] +
    agg.uncovered_factor * score[3] + agg.covered_factor * score[4]
=#

# Sums significance (1st), uncovered (3rd) and covered (4th) components of the score into a single score (`a`)
# if `-b/a` (b = 2nd component) is below ratio_threshold, keeps the score as is (`(a, b)`)
# otherwise, mixes in `b` into `a` and `a` into `b` (`(a + k sXs b, (1-k) b)`) with an increasing strength `k`,
# so that when `-b/a → ∞`, `(a#, b#) → (0, a/setXset_factor + b)`
# The transform keeps the aggregated score intact
struct MultiobjProblemSoftFold2d <: FitnessFolding{2}
    setXset_factor::Float64
    uncovered_factor::Float64
    covered_factor::Float64
    ratio_threshold::Float64    # at what set_score/setXset_score ratio the folding kicks in
    shape::Float64              # controls folding smoothness (smaller are smoother)

    MultiobjProblemSoftFold2d(params::CoverParams,
                              ratio_threshold::Float64 = 5.0,
                              shape::Float64=0.05) =
        new(params.setXset_factor, params.uncovered_factor, params.covered_factor,
            ratio_threshold, shape)
end

function (f::MultiobjProblemSoftFold2d)(fitness::RawScore)
    a = fitness[1] +
        f.uncovered_factor * fitness[3] +
        f.covered_factor * fitness[4]
    b = fitness[2]
    if f.ratio_threshold !== nothing
        t = b/max(1, -a)
        (t <= f.ratio_threshold) && return (a, b)

        ts = f.shape * (t - f.ratio_threshold)
        # 1/(1 + exp(1/ts - ts)) allows smooth transition from ts = 0 (f(0)=0) to f(∞)=1
        k = 1 / (t * f.setXset_factor) / (1 + exp(1/ts - ts))
        return (a + k * f.setXset_factor * b, (1.0 - k) * b)
    else
        return (a, b)
    end
end

(agg::MultiobjCoverProblemScoreAggregator{MultiobjProblemSoftFold2d})(score::NTuple{2, Float64}) =
    score[1] + agg.setXset_factor * score[2]

#=
# Sums uncovered (3rd) and covered (4th) components of the score into a single score
# that to every score component mixes in the other 2 components with the weight k
# The transform keeps the aggregated score intact
struct MultiobjProblemFold3d <: FitnessFolding{3}
    k::Float64
    k_sXs::Float64
    k_u::Float64
    setXset_factor::Float64
    uncovered_factor::Float64
    covered_factor::Float64

    MultiobjProblemFold3d(params::CoverParams, k::Float64=0.1) =
        new(params.setXset_factor > 0.0 ? k : 0.0,
            params.setXset_factor > 0.0 ? k/params.setXset_factor : 0.0,
            params.uncovered_factor > 0.0 ? k/params.uncovered_factor : 0.0,
            params.setXset_factor,
            params.uncovered_factor,
            params.uncovered_factor > 0.0 ? params.covered_factor/params.uncovered_factor : 1.0)
end

function (f::MultiobjProblemFold3d)(fitness::RawScore)
    a = fitness[1]
    b = fitness[2]
    c = fitness[3] + f.covered_factor * fitness[4]
    return (0.5 * f.k * (f.setXset_factor * b + f.uncovered_factor * c) + (1-f.k) * a,
            0.5 * f.k_sXs * (a + f.uncovered_factor * c) + (1-f.k) * b,
            0.5 * f.k_u * (a + f.setXset_factor * b) + (1-f.k) * c)
end

(agg::MultiobjCoverProblemScoreAggregator{MultiobjProblemFold3d})(score::NTuple{3, Float64}) =
    score[1] + agg.setXset_factor * score[2] + agg.uncovered_factor * score[3
=#

"""
Multi-objective optimal Enriched-Set Cover problem.
"""
struct MultiobjCoverProblem <: AbstractCoverProblem{RawScore}
    params::CoverParams

    var2set::Vector{Int}

    var_scores::Vector{Float64}
    varXvar_scores::Matrix{Float64}

    tileXvar::SparseMaskMatrix      # tile-in-mask X set(var)
    nmasked_tile::Vector{Int}       # sum of masked elements for each tile in all masks (els counted for each mask separately)
    nunmasked_tile::Vector{Int}     # number of unmasked elements for each tile in the union of all masks (each el counted once)

    tilepool::Vector{Vector{Float64}}   # pool of tile-sized vectors FIXME move to evaluator
    varpool::Vector{Vector{Float64}}    # pool of var-sized vectors FIXME move to evaluator

    function MultiobjCoverProblem(params::CoverParams,
                        var2set::AbstractVector{Int},
                        var_scores::AbstractVector{Float64},
                        varXvar_scores::AbstractMatrix{Float64},
                        tileXvar::SparseMaskMatrix,
                        nmasked_tile::AbstractVector{Int},
                        nunmasked_tile::AbstractVector{Int}
    )
        length(var2set) == length(var_scores) ==
        size(varXvar_scores, 1) == size(varXvar_scores, 2) ==
        size(tileXvar, 2) ||
            throw(ArgumentError("var2set, var_scores and varXvar_scores counts do not match"))
        size(tileXvar, 1) == length(nmasked_tile) == length(nunmasked_tile) ||
            throw(ArgumentError("tileXvar and nmasked_tile rows count does not match"))
        new(params, var2set,
            var_scores, varXvar_scores,
            tileXvar, nmasked_tile, nunmasked_tile,
            Vector{Vector{Float64}}(), Vector{Vector{Float64}}())
    end
end

ntiles(problem::MultiobjCoverProblem) = length(problem.nmasked_tile)

function score_scales(problem::MultiobjCoverProblem)
    w_min, w_max = extrema(problem.var_scores)
    sXs_min, sXs_max = length(problem.varXvar_scores) > 1 ? (Inf, -Inf) : (0.0, 0.0)
    return (max(w_max - w_min, 0.1),
            max(sXs_max - sXs_min, 0.1))
end

struct MultiobjOptimizerParams <: AbstractOptimizerParams{MultiobjCoverProblem}
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
        new(PopulationSize, WeightDigits, MaxSteps, MaxStepsWithoutProgress,
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
    tXv, nmasked_tile, nunmasked_tile = tilemaskXvar(mosaic)
    MultiobjCoverProblem(params,
                         v2set, v_scores, vXv_scores,
                         tXv, nmasked_tile, nunmasked_tile)
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
                problem.nmasked_tile, problem.nunmasked_tile)
end

# the tuple of:
# - total number of uncovered masked elements (in all masks overlapping with the cover)
# -`` total number of covered unmasked elements (in all masks overlapping with the cover)
function miscover_score(w::AbstractVector{Float64}, problem::MultiobjCoverProblem)
    __check_vars(w, problem)
    isempty(problem.nmasked_tile) && return (0.0, 0.0)
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
    res = sum(problem.nmasked_tile) - dot(problem.nmasked_tile, wtile),
          dot(problem.nunmasked_tile, wtile)
    push!(problem.tilepool, wtile) # release to the pool
    return res
end

"""
Unfolded multiobjective score (fitness) of the OESC coverage.

* `w` probabilities of the sets being covered
"""
function score(w::AbstractVector{Float64}, problem::MultiobjCoverProblem)
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
        c, d = miscover_score(w, problem)
    end
    return (a, b, c, d)
end

aggscore(s::RawScore, params::CoverParams) =
    s[1] +
    params.setXset_factor * s[2] +
    params.uncovered_factor * s[3] +
    params.covered_factor * s[4]

aggscore(w::AbstractVector{Float64}, problem::MultiobjCoverProblem) =
    aggscore(score(w, problem), problem.params)

const MultiobjCoverProblemFitnessScheme =
    ParetoFitnessScheme{2,Float64,true,MultiobjCoverProblemScoreAggregator{MultiobjProblemSoftFold2d}}

# wraps MultiobjCoverProblem as BlackBoxOptim OptimizationProblem
struct MultiobjCoverProblemBBOWrapper{SS <: RectSearchSpace} <:
        BlackBoxOptim.OptimizationProblem{MultiobjCoverProblemFitnessScheme}
    orig::MultiobjCoverProblem
    fitness_folding::MultiobjProblemSoftFold2d
    fitness_scheme::MultiobjCoverProblemFitnessScheme

    search_space::SS

    function MultiobjCoverProblemBBOWrapper(
        orig::MultiobjCoverProblem,
        fitness_folding::MultiobjProblemSoftFold2d,
        dimdigits::Union{Integer,Nothing}=2
    ) where FF
        fitness_scheme = MultiobjCoverProblemFitnessScheme(aggregator=MultiobjCoverProblemScoreAggregator{MultiobjProblemSoftFold2d}(orig.params))
        ss = RectSearchSpace(nvars(orig), (0.0, 1.0), dimdigits=dimdigits)
        new{typeof(ss)}(orig, fitness_folding, fitness_scheme, ss)
    end
end

Base.copy(problem::MultiobjCoverProblemBBOWrapper) =
    MultiobjCoverProblemBBOWrapper(problem.orig, problem.fitness_folding,
                                   dimdigits(search_space(problem), 1))
#=
BlackBoxOptim.show_fitness(io::IO, score::RawScore,
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
    @printf(io, "(sets+k·sets²=%.3f sets²+k·sets=%.3f)\n", score[1], score[2])

BlackBoxOptim.show_fitness(io::IO, score::IndexedTupleFitness, problem::MultiobjCoverProblemBBOWrapper) =
    BlackBoxOptim.show_fitness(io, score.orig, problem)

BlackBoxOptim.fitness(x::BlackBoxOptim.AbstractIndividual,
                      p::MultiobjCoverProblemBBOWrapper) = p.fitness_folding(score(x, p.orig))

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
Result of `optimize(MultiobjCoverProblem)`.
"""
struct MultiobjCoverProblemResult
    var2set::Vector{Int}        # indices of sets in the original SetMosaic
    varscores::Vector{Float64}  # scores of the sets

    varweights::Matrix{Float64} # weights of the vars for each solution on Pareto front
    raw_scores::Vector{RawScore}# raw scores for each solution
    folded_scores::Vector{FoldedScore}
    agg_scores::Vector{Float64}
    best_ix::Int                # index of the best solution (agg_score minimum)
end

function MultiobjCoverProblemResult(problem::MultiobjCoverProblem,
                                    optres::BlackBoxOptim.OptimizationResults)
    # initialize solutions
    wmtx = similar(best_candidate(optres), nvars(problem),
                   length(pareto_frontier(optres)))
    folded_scores = [archived_fitness(sol).orig for sol in pareto_frontier(optres)]
    front_perm = sortperm(folded_scores) # sort lexicographically
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
    # calculate raw scores for the corrected weights
    raw_scores = dropdims(mapslices(w -> score(w, problem), wmtx; dims=1), dims=1)
    agg_scores = aggscore.(raw_scores, Ref(problem.params))

    return MultiobjCoverProblemResult(problem.var2set, problem.var_scores,
                                      wmtx, raw_scores, folded_scores, agg_scores,
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
    min_score = aggscore(res.raw_scores[best_ix], params)
    @inbounds for i in 2:nsolutions(res)
        sol_score = aggscore(res.raw_scores[i], params)
        if sol_score < min_score
            best_ix = i
            min_score = sol_score
        end
    end
    return best_ix
end

best_varweights(result::MultiobjCoverProblemResult) =
    @inbounds(view(result.varweights, :, best_index(result)))

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

function optimize(problem::MultiobjCoverProblem,
                  opt_params::MultiobjOptimizerParams = MultiobjOptimizerParams())
    if nvars(problem) == 0
        return MultiobjCoverProblemResult(
            Vector{Int}(), Vector{Float64}(), Matrix{Float64}(undef, 0, 0),
            Vector{RawScore}(), Vector{FoldedScore}(), Vector{Float64}(), 0)
    end

    fitfolding = MultiobjProblemSoftFold2d(problem.params)
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
    return MultiobjCoverProblemResult(problem, bbores)
end
