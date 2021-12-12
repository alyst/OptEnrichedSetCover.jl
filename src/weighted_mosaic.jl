"""
[`SetMosaic`](@ref) with the weights assigned to its sets.

## Type parameters
* `T`: the type of elements
* `S`: the type of set ids
* `E`: the type of experiment ids
* `W`: type of the weight
"""
abstract type AbstractWeightedSetMosaic{T,S,E,W} end

"""
    originalmosaic(mosaic::AbstractWeightedSetMosaic) -> SetMosaic

Get the original [`SetMosaic`](@ref).
"""
originalmosaic(mosaic::AbstractWeightedSetMosaic) = mosaic.original

Base.eltype(::Type{<:AbstractWeightedSetMosaic{T}}) where T = T
Base.eltype(mosaic::AbstractWeightedSetMosaic) = eltype(typeof(mosaic))

setidtype(::Type{<:AbstractWeightedSetMosaic{<:Any, S}}) where S = S
setidtype(mosaic::AbstractWeightedSetMosaic) = setidtype(typeof(mosaic))

experimentidtype(::Type{<:AbstractWeightedSetMosaic{<:Any, <:Any, E}}) where E = E
experimentidtype(mosaic::AbstractWeightedSetMosaic) = experimentidtype(typeof(mosaic))

weighttype(::Type{<:AbstractWeightedSetMosaic{<:Any, <:Any, <:Any, W}}) where W = W
weighttype(mosaic::AbstractWeightedSetMosaic) = weighttype(typeof(mosaic))

"""
[`SetMosaic`](@ref) with the weights for the sets from multiple experiments on top.

## Type parameters
* `T`: the type of elements
* `S`: the type of set ids
* `E`: the type of experiment ids
"""
mutable struct WeightedSetMosaic{T,S,E,W} <: AbstractWeightedSetMosaic{T,S,E,W}
    original::SetMosaic{T,S}        # original mosaic

    setXexp_weights::Matrix{W}      # setXexperiment weights matrix
    max_weight::Union{W, Nothing}

    loc2glob_setix::Vector{Int}     # local set index to global set index in original mosaic
    glob2loc_setix::SparseVector{Int}   # global set index (in orig. mosaic) to local set index

    ix2experiment::Vector{E}        # experiment index to the ID of the experiment
    experiment2ix::Dict{E, Int}     # experiment ID to index
end

function WeightedSetMosaic(mosaic::SetMosaic{T, S},
                           setXexp_weights::AbstractMatrix{W},
                           max_weight::Union{Number, Nothing},
                           set_indices::Vector{Int},
                           experiment_ids::Union{AbstractVector{E}, Base.KeySet{E}}
                          ) where {T,S,W,E}
    size(setXexp_weights, 1) == length(set_indices) ||
        throw(ArgumentError("Set weights rows ($(size(setXexp_weights, 1))) should match the number of sets ($(length(set_indices)))"))
    nexps = size(setXexp_weights, 2)
    (experiment_ids === nothing) || (nexps == length(experiment_ids)) ||
        throw(ArgumentError("Number of masks ($nexps) should match the number of mask IDs ($(length(experiment_ids)))"))
    WeightedSetMosaic(mosaic, setXexp_weights,
        !isnothing(max_weight) ? convert(nonmissingtype(W), max_weight) : nothing,
        set_indices, SparseVector(nsets(mosaic), set_indices, collect(eachindex(set_indices))),
        experiment_ids isa AbstractVector ?
            convert(Vector, experiment_ids) :
            collect(experiment_ids),
        # FIXME use vector container if experiment_ids === nothing or if ids are integers
        Dict(id => ix for (ix, id) in enumerate(experiment_ids)))
end

"""
    assignweights(mosaic::SetMosaic, elmasks::AbstractMatrix{Bool};
                 [experiment_ids::Union{AbstractVector, AbstractSet, Nothing} = nothing],
                 [max_setsize=nothing],
                 [max_weight], [max_min_weight]) -> WeightedSetMosaic

Construct [`WeightedSetMosaic`](@ref) from the [`SetMosaic`](@ref)
and the external `weights`.

## Arguments
* `max_setsize` (optional): ignore the annotation sets bigger than the specified size
* `max_weight` (optional): the maximal weight of the set to include in the mosaic
* `max_min_weight` (optional): the maximual weight of the set *in all experiments* to include the set into mosaic
"""
function assignweights(mosaic::SetMosaic, weights::AbstractMatrix;
                       set_ids::Union{AbstractVector, Nothing} = nothing,
                       experiment_ids::Union{AbstractVector, Base.KeySet, Nothing} = nothing,
                       max_setsize::Union{Integer, Nothing} = nothing,
                       max_weight::Union{Number, Nothing} = nothing,
                       max_min_weight::Union{Number, Nothing} = nothing
)
    __set_ids = isnothing(set_ids) ? mosaic.ix2set : set_ids
    size(weights, 1) == length(__set_ids) ||
        throw(ArgumentError("Weight rows ($(size(weights, 1))) should match the number of set IDs ($(length(__set_ids)))"))
    isnothing(experiment_ids) || size(weights, 2) == length(experiment_ids) ||
        throw(ArgumentError("Weight columns ($(size(weights, 2))) should match the number of experiment IDs ($(length(experiment_ids)))"))

    set2row = Dict(s => i for (i, s) in enumerate(__set_ids))
    # get the sets that pass the weight filters
    setixs = Vector{Int}()
    for (setix, s) in enumerate(mosaic.ix2set)
        rowix = get(set2row, s, 0)
        (rowix == 0) && continue # no weights provided
        if !isnothing(max_setsize) && setsize(mosaic, setix) > max_setsize
            continue # too large
        end
        rowweights = view(weights, rowix, :)
        all(ismissing, rowweights) && continue # no weights provided
        !isnothing(max_min_weight) &&
            (minimum(skipmissing(rowweights)) > max_min_weight) && continue # too high weights
        push!(setixs, setix)
    end
    # calculate weight type
    W = nonmissingtype(eltype(weights))
    # check if it's necessary to store missing values
    if any(w -> ismissing(w) || !isnothing(max_weight) && (w > max_weight), vec(weights))
        W = Union{W, Missing}
    end
    # create the setXexp_weights from weights matrix
    setXexp_weights = Matrix{W}(undef, length(setixs), size(weights, 2))
    for (i, setix) in enumerate(setixs)
        rowix = set2row[mosaic.ix2set[setix]]
        setXexp_weights[i, :] .= view(weights, rowix, :)
    end

    return WeightedSetMosaic(mosaic, setXexp_weights, max_weight, setixs, experiment_ids)
end

function _assignweights(mosaic::SetMosaic{<:Any, S},
                        exps_weights::Union{AbstractVector, Base.ValueIterator}, #= iterable with eltype()==Dict{T, W} =#
                        exps_ids::Union{AbstractVector, Base.KeySet};
                        kwargs...) where S
    @assert eltype(exps_weights) <: AbstractDict{S}
    @assert length(exps_weights) == length(exps_ids)
    allsets = Set{S}()
    for exp_weights in exps_weights
        for s in keys(exp_weights)
            push!(allsets, s)
        end
    end
    # indices of all sets that exist in the mosaic
    set_ixs = get.(Ref(mosaic.set2ix), allsets, 0)
    nmissingsets = isempty(set_ixs) ? 0 : sum(==(0), set_ixs)
    (nmissingsets > 0) && @warn "$nmissingsets set(s) not found in set mosaic"

    # drop sets that are not found in the mosaic
    sort!(filter!(!=(0), set_ixs))
    set2ix = Dict(mosaic.ix2set[ix] => i for (i, ix) in enumerate(set_ixs))
    # convert exps_weights to a matrix
    W = nonmissingtype(valtype(eltype(exps_weights)))
    weights_mtx = fill!(Matrix{Union{W, Missing}}(undef, length(set_ixs), length(exps_weights)), missing)
    for (expix, exp_weights) in enumerate(exps_weights)
        exp_weights_v = view(weights_mtx, :, expix)
        for (s, w) in pairs(exp_weights)
            setix = get(set2ix, s, 0)
            setix == 0 && continue # skip weight of missing set
            exp_weights_v[setix] = w
        end
    end

    assignweights(mosaic, weights_mtx,
                  set_ids = mosaic.ix2set[set_ixs],
                  experiment_ids = exps_ids; kwargs...)
end

assignweights(mosaic::SetMosaic,
              weights::AbstractVector; #= iterable with eltype()==Dict{T, W} =#
              kwargs...) =
    _assignweights(mosaic, weights, eachindex(weights); kwargs...)

assignweights(mosaic::SetMosaic,
              weights::AbstractDict; #= iterable with eltype()==Dict{T, W} =#
              kwargs...) =
    _assignweights(mosaic, values(weights), keys(weights); kwargs...)

nelements(mosaic::WeightedSetMosaic) = nelements(originalmosaic(mosaic))
nsets(mosaic::WeightedSetMosaic) = length(mosaic.loc2glob_setix) # only sets overlapping with masks
nexperiments(mosaic::WeightedSetMosaic) = length(mosaic.ix2experiment)

function setweight(mosaic::WeightedSetMosaic, set::Any, experiment::Any;
                   filter::Bool = true)
    expix = mosaic.experiment2ix[experiment]
    glob_setix = mosaic.original.set2ix[set]
    setix = mosaic.glob2loc_setix[glob_setix]
    return setix != 0 ? _setweight(mosaic, setix, expix, filter=filter) : missing
end

function _setweight(mosaic::WeightedSetMosaic, setix::Integer, expix::Integer;
                    filter::Bool = true)
    w = mosaic.setXexp_weights[setix, expix]
    ifelse(!filter || isnothing(mosaic.max_weight) ||
           (!ismissing(w) && (w <= mosaic.max_weight)),
           w, missing)
end

# copy everything, except the original mosaic (leave the reference to the same object)
Base.copy(mosaic::WeightedSetMosaic) =
    WeightedSetMosaic(mosaic.original, copy(mosaic.setXexp_weights), mosaic.max_weight,
                      copy(mosaic.loc2glob_setix), copy(mosaic.glob2loc_setix),
                      copy(mosaic.ix2experiment), copy(mosaic.experiment2ix))
