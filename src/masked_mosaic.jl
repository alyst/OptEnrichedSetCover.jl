# the number of masked elements in each set of mosaic
function nmasked_perset_mask!(nmasked::AbstractVector{Int}, mosaic::SetMosaic, elmask::AbstractVector{Bool},
                              set2index::Union{Nothing,Dict{Int,Int}} = nothing)
    @inbounds for elix in eachindex(elmask)
        elmask[elix] || continue
        for setix in view(mosaic.setXelm, :, elix)
            if set2index === nothing
                nmasked[setix] += 1
            else
                ix = get(set2index, setix, 0)
                if ix > 0
                    nmasked[ix] += 1
                end
            end
        end
    end
    return nmasked
end

# the number of masked elements in each set of mosaic
function nmasked_perset(mosaic::SetMosaic, elmasks::AbstractVector,
                        set2index::Union{Nothing,Dict{Int,Int}} = nothing)
    res = zeros(Int, set2index === nothing ? nsets(mosaic) : length(set2index), length(elmasks))
    @inbounds for (expix, elmask) in enumerate(elmasks)
        @assert length(elmask) == nelements(mosaic)
        nmasked_perset_mask!(view(res, :, expix), mosaic, elmask, set2index)
    end
    return res
end

# the number of masked elements in each set of mosaic
function nmasked_perset(mosaic::SetMosaic, elmasks::AbstractMatrix{Bool},
                        set2index::Union{Nothing,Dict{Int,Int}} = nothing)
    @assert size(elmasks, 1) == nelements(mosaic)
    res = zeros(Int, set2index === nothing ? nsets(mosaic) : length(set2index), size(elmasks, 2))
    @inbounds for expix in 1:size(elmasks, 2)
        nmasked_perset_mask!(view(res, :, expix), mosaic, view(elmasks, :, expix), set2index)
    end
    return res
end

# info for the set × mask overlap for the MaskedSetMosaic
struct MaskOverlap
    expix::Int      # index of the experiment
    nmasked::Int    # number of masked elements
    nunmasked::Int  # number of unmasked elements
end

"""
[`SetMosaic`](@ref) with the elements masks (selections) on top.
Sets that are not overlapping with the masks are excluded(skipped) from `MaskedSetMosaic`.
Optionally, the filtering can include testing for the minimal overlap significance P-value.

The tiles of non-overlapped sets are removed, the tiles that have identical membership
for all the masked sets are squashed into a single tile.

## Type parameters
* `T`: the type of elements
* `S`: the type of set ids
* `E`: the type of experiment ids
"""
mutable struct MaskedSetMosaic{T, S, E}
    original::SetMosaic{T, S}       # original mosaic
    elmasks::BitMatrix              # elements×masks
    elunionmask::BitVector          # all masked elements
    total_masked::Vector{Int}       # total number of masked elements for each mask

    # info for all non-empty overlaps with experiment masks
    # experiments in MaskOverlap vector are ordered by experiment index
    set2experiments::Dict{Int, Vector{MaskOverlap}}

    ix2experiment::Vector{E}        # experiment index to the ID of the mask
    experiment2ix::Dict{E, Int}     # experiment Id to index
end

function MaskedSetMosaic(mosaic::SetMosaic{T, S}, elmasks::AbstractMatrix{Bool},
                         set2experiments::Dict{Int, Vector{MaskOverlap}},
                         experiment_ids::Union{AbstractVector{E}, Nothing} = nothing
                        ) where {T, S, E}
    size(elmasks, 1) == nelements(mosaic) ||
        throw(ArgumentError("Elements mask length ($(size(elmasks, 1))) should match the number of elements ($(nelements(mosaic)))"))
    nmasks = size(elmasks, 2)
    (experiment_ids === nothing) || (nmasks == length(experiment_ids)) ||
        throw(ArgumentError("Number of masks ($nmasks) should match the number of mask IDs ($(length(masks_ids)))"))
    bit_elmasks = convert(BitMatrix, elmasks)
    bit_elunionmask = dropdims(any(bit_elmasks, dims=2), dims=2)
    @assert length(bit_elunionmask) == size(bit_elmasks, 1)
    MaskedSetMosaic(mosaic, bit_elmasks, bit_elunionmask,
                    dropdims(sum(elmasks, dims=1), dims=1), set2experiments,
                    collect(experiment_ids !== nothing ? experiment_ids : 1:nmasks),
        # FIXME use vector container if experiment_ids === nothing or if ids are integers
        Dict(id => ix for (ix, id) in enumerate(experiment_ids !== nothing ? experiment_ids : 1:nmasks)))
end

"""
    mask(mosaic::SetMosaic, elmasks::AbstractMatrix{Bool};
         [experiment_ids::Union{AbstractVector, AbstractSet, Nothing} = nothing],
         [min_nmasked=1], [max_setsize=nothing],
         [max_overlap_logpvalue=0.0]) -> MaskedSetMosaic

Construct [`MaskedSetMosaic`](@ref) from the [`SetMosaic`](@ref) and the collection of element masks.

## Arguments
* `min_nmasked`: the minimal number of masked elements in a set to include in the mosaic
* `max_setsize` (optional): ignore the annotation sets bigger than the specified size
* `max_overlap_logpvalue`: the threshold of Fisher's Exact Test log P-value of the overlap
   between the set and the mask for the inclusion of the set into the mosaic.
   *0* accepts all sets.
"""
function mask(mosaic::SetMosaic, elmasks::AbstractMatrix{Bool};
              experiment_ids::Union{AbstractVector, Nothing} = nothing,
              min_nmasked::Integer=1, max_setsize::Union{Integer, Nothing} = nothing,
              max_overlap_logpvalue::Number=0.0, # 0.0 would accept any overlap (as log(Fisher Exact Test P-value))
              max_min_overlap_logpvalue::Number=max_overlap_logpvalue
)
    size(elmasks, 1) == nelements(mosaic) ||
        throw(ArgumentError("Elements mask rows ($(size(elmasks, 1))) should match the number of elements ($(nelements(mosaic)))"))

    # calculate overlap sizes for each set and each mask
    nmasked_orgsets = nmasked_perset(mosaic, elmasks)

    # get the sets that overlap with the masked elements and with at least max_overlap_logpvalue significance
    set2experiments = Dict{Int, Vector{MaskOverlap}}()
    emptyoverlap = Vector{MaskOverlap}()
    ntotal = nelements(mosaic)
    min_overlap_logpvalue = fill(0.0, size(nmasked_orgsets, 1))
    @inbounds for expix in axes(nmasked_orgsets, 2)
        ntotal_masked = sum(view(elmasks, :, expix))
        orgsets_mask = view(nmasked_orgsets, :, expix)
        for (setix, nmasked) in enumerate(orgsets_mask)
            (nmasked < min_nmasked) && continue # skip small overlap
            @inbounds nset = setsize(mosaic, setix)
            (max_setsize !== nothing) && (nset > max_setsize) && continue # skip very big and generic sets
            overlap_pvalue = logpvalue(nmasked, nset, ntotal_masked, ntotal)
            (overlap_pvalue > max_overlap_logpvalue) && continue # skip non-signif overlaps
            (overlap_pvalue < min_overlap_logpvalue[setix]) && (min_overlap_logpvalue[setix] = overlap_pvalue)
            setmasks = get!(() -> Vector{MaskOverlap}(), set2experiments, setix)
            push!(setmasks, MaskOverlap(expix, nmasked, setsize(mosaic, setix) - nmasked))
        end
    end
    filter!(kv -> min_overlap_logpvalue[kv[1]] <= max_min_overlap_logpvalue, pairs(set2experiments))

    return MaskedSetMosaic(mosaic, elmasks, set2experiments, experiment_ids)
end

function mask(mosaic::SetMosaic{T}, elmasks::AbstractVector #= iterable with eltype()==Set{T} =#;
              kwargs...) where T
    @assert eltype(elmasks) === Set{T}
    mask(mosaic, Bool[in(e, elmask::Set{T}) for e in mosaic.ix2elm, elmask in elmasks]; kwargs...)
end

function mask(mosaic::SetMosaic{T}, elmasks::AbstractDict #= iterable with eltype()==Set{T} =#;
              kwargs...) where T
    @assert valtype(elmasks) === Set{T}
    mask(mosaic, Bool[in(e, elmask::Set{T}) for e in mosaic.ix2elm, elmask in values(elmasks)];
         experiment_ids = collect(keys(elmasks)), kwargs...)
end

"""
    originalmosaic(mosaic::MaskedSetMosaic) -> SetMosaic

Get the original [`SetMosaic`](@ref).
"""
originalmosaic(mosaic::MaskedSetMosaic) = mosaic.original

nelements(mosaic::MaskedSetMosaic) = nelements(mosaic.original)
nsets(mosaic::MaskedSetMosaic) = length(mosaic.set2experiments) # only sets overlapping with masks
nsets(mosaic::MaskedSetMosaic, expix::Int) = # number of sets overlapping with given mask
    mapreduce(olaps -> any(olap -> olap.mask == expix), sum, values(mosaic.set2experiments), init=0)
nmasks(mosaic::MaskedSetMosaic) = length(mosaic.total_masked)
_nmasked(mosaic::MaskedSetMosaic, expix::Int) = mosaic.total_masked[expix]
_nunmasked(mosaic::MaskedSetMosaic, expix::Int) = nelements(mosaic) - mosaic.total_masked[expix]

nmasked(mosaic::MaskedSetMosaic, mask::Any) = _nmasked(mosaic, mosaic.experiment2ix[mask])
nunmasked(mosaic::MaskedSetMosaic, mask::Any) = _nunmasked(mosaic, mosaic.experiment2ix[mask])

# FIXME swap exp and set args order
function overlap(mosaic::MaskedSetMosaic, exp::Any, set::Any)
    expix = mosaic.experiment2ix[exp]
    setix = mosaic.original.set2ix[set]
    olaps = get(mosaic.set2experiments, setix, nothing)
    isnothing(olaps) && return nothing
    olapix = searchsortedfirst(olaps, MaskOverlap(expix, -1, -1), by=x -> x.mask)
    ((olapix > length(olaps)) || (olaps[olapix].mask != expix)) && return nothing
    return olaps[olapix]
end

# FIXME swap mask and set args order
function logpvalue(mosaic::MaskedSetMosaic, exp::Any, set::Any)
    olap = overlap(mosaic, exp, set)
    isnothing(olap) && return missing
    return logpvalue(olap.nmasked, olap.nmasked + olap.nunmasked,
                     nmasked(mosaic, exp), nelements(mosaic))
end

# copy everything, except the original mosaic (leave the reference to the same object)
Base.copy(mosaic::MaskedSetMosaic) =
    MaskedSetMosaic(mosaic.original, copy(mosaic.elmasks), copy(mosaic.elunionmask),
                    copy(mosaic.total_masked),
                    deepcopy(mosaic.set2experiments),
                    copy(mosaic.ix2experiment), copy(mosaic.experiment2ix))
