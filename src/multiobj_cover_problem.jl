"""
Multi-objective optimal Enriched-Set Cover problem.
"""
struct MultiobjectiveCoverProblem <: AbstractCoverProblem{NTuple{3, Float64}}
    params::CoverParams

    var2mset::Vector{Int}
    var_ranges::Vector{UnitRange}           # variables grouped by mask
    var_scores::Vector{Float64}
    varXvar_scores::Matrix{Matrix{Float64}}

    function MultiobjectiveCoverProblem(params::CoverParams,
                        var2mset::AbstractVector{Int},
                        var_ranges::AbstractVector{UnitRange},
                        var_scores::AbstractVector{Float64},
                        varXvar_scores::AbstractMatrix{Matrix{Float64}})
        length(var_ranges) ==
        size(varXvar_scores, 1) == size(varXvar_scores, 2) ||
            throw(ArgumentError("var_scores and varXvar_scores block counts do not match"))
        if !isempty(var2mset)
            length(var_scores) == length(var2mset) == sum(length, var_ranges) ||
                throw(ArgumentError("var_scores and var2mset sizes do not match"))
            length.(var_ranges) ==
            size.(view(varXvar_scores, :, 1), 1) == size.(view(varXvar_scores, 1, :), 2) ||
                throw(ArgumentError("var_scores and varXvar_scores block sizes do not match"))
        else
            isempty(var_scores) || throw(ArgumentError("var_scores and var2mset sizes do not match"))
        end
        new(params, var2mset, var_ranges,
            var_scores, varXvar_scores)
    end
end

nmasks(problem::MultiobjectiveCoverProblem) = length(problem.var_ranges)

function score_scales(problem::MultiobjectiveCoverProblem)
    w_min, w_max = extrema(problem.var_scores)
    sXs_min, sXs_max = Inf, -Inf
    mXm_min, mXm_max = length(problem.varXvar_scores) > 1 ? (Inf, -Inf) : (0.0, 0.0)
    for i in 1:size(problem.varXvar_scores, 1)
        isempty(problem.var_ranges[i]) && continue # skip empty (no extrema)
        sXs_min_i, sXs_max_i = extrema(problem.varXvar_scores[i, i])
        (sXs_min > sXs_min_i) && (sXs_min = sXs_min_i)
        (sXs_max < sXs_max_i) && (sXs_max = sXs_max_i)
        for j in 1:size(problem.varXvar_scores, 2)
            ((i == j) || isempty(problem.var_ranges[j])) && continue
            mXm_min_ij, mXm_max_ij = extrema(problem.varXvar_scores[i, j])
            (mXm_min > mXm_min_ij) && (mXm_min = mXm_min_ij)
            (mXm_max < mXm_max_ij) && (mXm_max = mXm_max_ij)
        end
    end
    return (max(w_max - w_min, 0.1),
            max(sXs_max - sXs_min, 0.1),
            max(mXm_max - mXm_min, 0.1))
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
        eps_scale = get(opt_params.borg_params, :ϵ_scale, [0.2, 0.1, 0.1])
        eps = eps_scale .* [score_scales(problem)...]
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

function MultiobjectiveCoverProblem(mosaic::MaskedSetMosaic, params::CoverParams = CoverParams())
    var2mset = collect(1:nsets(mosaic))
    # calculate variable ranges for each mask
    var_ranges = Vector{UnitRange}()
    maskset1 = nextset = 1
    while nextset <= length(mosaic.maskedsets)
        nextset += 1
        curmask = mosaic.maskedsets[maskset1].mask
        nextmask = nextset > length(mosaic.maskedsets) ?
                   nmasks(mosaic)+1 : mosaic.maskedsets[nextset].mask
        if curmask != nextmask
            @assert curmask < nextmask
            # insert empty masks before curmask
            while length(var_ranges) < curmask-1
                push!(var_ranges, maskset1:(maskset1-1))
                @assert nsets(mosaic, length(var_ranges)) == 0
            end
            push!(var_ranges, maskset1:(nextset-1))
            maskset1 = nextset
        end
    end
    # fill with empty masks to the end
    while length(var_ranges) < nmasks(mosaic)
        push!(var_ranges, maskset1:(maskset1-1))
        @assert nsets(mosaic, length(var_ranges)) == 0
    end
    var_scores = msetscore_detached.(mosaic.maskedsets, mosaic, params) .- log(params.sel_prob)

    # prepare varXvar scores block-matrix
    varXvar_scores = Matrix{typeof(mosaic.original.setXset_scores)}(length(var_ranges), length(var_ranges))
    min_sXs = Inf # minimum finite varXvar_scores element
    @inbounds for (ii, irange) in enumerate(var_ranges)
        for jj in ii:length(var_ranges)
            jrange = var_ranges[jj]
            varXvar_scores[ii, jj] = mXm_mtx = similar(mosaic.original.setXset_scores, length(irange), length(jrange))
            for (i, iset) in enumerate(view(mosaic.maskedsets, irange))
                for (j, jset) in enumerate(view(mosaic.maskedsets, jrange))
                    sXs = mosaic.original.setXset_scores[iset.set, jset.set]
                    mXm_mtx[i, j] = varXvar_score(sXs, iset, jset, params, false)
                    if iset.set != jset.set && isfinite(sXs) && sXs < min_sXs
                        min_sXs = sXs
                    end
                end
            end
            (jj > ii) && (varXvar_scores[jj, ii] = transpose(mXm_mtx)) # the block-matrix is symmetric
        end
    end
    # replace infinite varXvar score with the minimal finite setXset score
    sXs_min = varXvar_score(1.25*min_sXs, MaskedSet(1, 1, 0, 0), MaskedSet(2, 1, 0, 0), params, true)
    mXm_min = varXvar_score(1.25*min_sXs, MaskedSet(1, 1, 0, 0), MaskedSet(2, 2, 0, 0), params, true)
    @inbounds for (ii, mXm_mtx) in enumerate(varXvar_scores)
        m1, m2 = ind2sub(size(varXvar_scores), ii)
        for i in eachindex(mXm_mtx)
            if !isfinite(mXm_mtx[i])
                v1, v2 = ind2sub(size(mXm_mtx), i)
                warn("var[$(var_ranges[m1][v1])]×var[$(var_ranges[m2][v2])] score is $(mXm_mtx[i])")
                if mXm_mtx[i] < 0.0
                    mXm_mtx[i] = m1 == m2 ? sXs_min : mXm_min
                end
            end
        end
    end
    MultiobjectiveCoverProblem(params, var2mset,
                               var_ranges, var_scores, varXvar_scores)
end

function exclude_vars(problem::MultiobjectiveCoverProblem,
                      vars::AbstractVector{Int};
                      penalize_overlaps::Bool = true)
    # constuct new var2mset and var_ranges
    newvar_ranges = similar(problem.var_ranges)
    var2mset = similar(problem.var2mset)
    varmask = fill(true, nvars(problem))
    varmask[vars] = false
    var2mset = problem.var2mset[varmask]
    last_mask_newvar = 0
    newvar_ranges = similar(problem.var_ranges)
    for (i, var_range) in enumerate(problem.var_ranges)
        range_len = sum(view(varmask, var_range))
        newvar_ranges[i] = (last_mask_newvar+1):(last_mask_newvar+range_len)
        last_mask_newvar += range_len
    end
    var_scores = copy(problem.var_scores)
    if penalize_overlaps
        penalty_weights = fill(0.0, nvars(problem))
        penalty_weights[vars] = 1.0
        wXw = fill!(similar(penalty_weights), 0.0) # intermask overlap penalties
        for i in 1:size(problem.varXvar_scores, 1)
            iseg = problem.var_ranges[i]
            wi = view(penalty_weights, iseg)
            tmpi = problem.varXvar_scores[i, i] * wi
            tmpi .*= wi
            # penalize within-mask overlaps
            var_scores_i = view(var_scores, iseg)
            var_scores_i .-= problem.params.setXset_factor .* tmpi
            # update intermask overlaps
            for j in 1:size(problem.varXvar_scores, 2)
                (i == j) && continue
                jseg = problem.var_ranges[j]
                view(wXw, jseg) .+= (problem.varXvar_scores[i, j] * view(penalty_weights, jseg)) .* wi
            end
        end
        # penalize intermask overlaps
        var_scores .-= (problem.params.setXset_factor*problem.params.maskXmask_factor) .* wXw
    end
    var_scores = var_scores[varmask] # filter out unused vars
    # filter out vars from the matrix
    varXvar_scores = similar(problem.varXvar_scores)
    for j in 1:size(problem.varXvar_scores, 2)
        jseg = problem.var_ranges[j]
        jmask = view(varmask, jseg)
        for i in 1:size(problem.varXvar_scores, 1)
            iseg = problem.var_ranges[i]
            imask = view(varmask, iseg)
            varXvar_scores[i, j] = problem.varXvar_scores[i, j][imask, jmask]
        end
    end
    return MultiobjectiveCoverProblem(problem.params,
                var2mset, newvar_ranges,
                var_scores, varXvar_scores)
end

"""
Score (probability) of the OESC coverage.

* `w` probabilities of the sets being covered
"""
function score(problem::MultiobjectiveCoverProblem, w::AbstractVector{Float64})
    if problem.params.setXset_factor == 0.0
        # skip set interactions as it would be aggregated to zero
        return (dot(problem.var_scores, w), 0.0, 0.0)
    end
    maskXmask = 0.0
    setXset = 0.0
    ww = fill!(similar(w), 0.0)
    for i in 1:size(problem.varXvar_scores, 1)
        iseg = problem.var_ranges[i]
        wi = view(w, iseg)
        wwi = fill(0.0, length(iseg))
        tmpi = similar(wwi)
        A_mul_B!(tmpi, problem.varXvar_scores[i, i], wi)
        setXset += dot(tmpi, wi)
        for j in 1:size(problem.varXvar_scores, 2)
            (i == j) && continue
            wj = view(w, problem.var_ranges[j])
            A_mul_B!(tmpi, problem.varXvar_scores[i, j], wj)
            wwi .= min.(wwi, tmpi)
        end
        maskXmask += dot(wwi, wi)
    end
    return (dot(problem.var_scores, w), -setXset, -maskXmask)
end

struct MultiobjectiveCoverProblemBBOWrapper <: BlackBoxOptim.OptimizationProblem{ParetoFitnessScheme{3}}
    orig::MultiobjectiveCoverProblem
    search_space::SearchSpace

    function MultiobjectiveCoverProblemBBOWrapper(orig::MultiobjectiveCoverProblem)
        new(orig, symmetric_search_space(nvars(orig), (0.0, 1.0)))
    end
end

struct MultiobjectiveCoverProblemScoreAggregator
    k_sXs::Float64
    k_mXm::Float64

    MultiobjectiveCoverProblemScoreAggregator(params::CoverParams) =
        new(params.setXset_factor, params.maskXmask_factor * params.setXset_factor)
end

(agg::MultiobjectiveCoverProblemScoreAggregator)(score::NTuple{3, Float64}) =
    score[1] + agg.k_sXs * score[2] + agg.k_mXm * score[3]

aggscore(problem::MultiobjectiveCoverProblem, w::AbstractVector{Float64}) =
    MultiobjectiveCoverProblemScoreAggregator(problem.params)(score(problem, w))

BlackBoxOptim.name(::MultiobjectiveCoverProblemBBOWrapper) = "MultiObjectiveOptimalSetCoverProblem"

BlackBoxOptim.numdims(p::MultiobjectiveCoverProblemBBOWrapper) = BlackBoxOptim.numdims(p.search_space)
BlackBoxOptim.search_space(p::MultiobjectiveCoverProblemBBOWrapper) = p.search_space
BlackBoxOptim.fitness_scheme(p::MultiobjectiveCoverProblemBBOWrapper) =
    ParetoFitnessScheme{3}(is_minimizing=true,
                           aggregator=MultiobjectiveCoverProblemScoreAggregator(p.orig.params))

Base.copy(problem::MultiobjectiveCoverProblemBBOWrapper) = MultiobjectiveCoverProblemBBOWrapper(problem.orig)

BlackBoxOptim.show_fitness(io::IO, score::NTuple{3,Float64}, problem::MultiobjectiveCoverProblemBBOWrapper) =
   @printf(io, "  sets=%.3f set×set=%.3f mask×mask=%.3f\n", score...)

BlackBoxOptim.show_fitness(io::IO, score::IndexedTupleFitness{2,Float64}, problem::MultiobjectiveCoverProblemBBOWrapper) =
    BlackBoxOptim.show_fitness(io, score.orig, problem)

BlackBoxOptim.fitness(x, p::MultiobjectiveCoverProblemBBOWrapper) = score(p.orig, x)

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
        return CoverProblemResult(Vector{Int}(), Vector{Float64}(), Vector{Float64}(),
                                  (0.0, 0.0, 0.0), 0.0, nothing)
    end

    bbowrapper = MultiobjectiveCoverProblemBBOWrapper(problem)
    popmatrix = BlackBoxOptim.rand_individuals_lhs(search_space(bbowrapper), opt_params.pop_size)
    # two extreme solutions and one totally neutral
    size(popmatrix, 2) > 0 && (popmatrix[:, 1] = 0.0)
    size(popmatrix, 2) > 1 && (popmatrix[:, 2] = 1.0)
    size(popmatrix, 2) > 2 && (popmatrix[:, 3] = 0.5)
    population = FitPopulation(popmatrix, nafitness(IndexedTupleFitness{3,Float64}), ntransient=1)

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

    return CoverProblemResult(problem.var2mset, w, problem.var_scores .* w, s,
                              BlackBoxOptim.aggregate(s, fitness_scheme(bbores)),
                              fitness_frontier)
end
