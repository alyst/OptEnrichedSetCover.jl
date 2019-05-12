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
    @inbounds for (maskix, elmask) in enumerate(elmasks)
        @assert length(elmask) == nelements(mosaic)
        nmasked_perset_mask!(view(res, :, maskix), mosaic, elmask, set2index)
    end
    return res
end

# the number of masked elements in each set of mosaic
function nmasked_perset(mosaic::SetMosaic, elmasks::AbstractMatrix{Bool},
                        set2index::Union{Nothing,Dict{Int,Int}} = nothing)
    @assert size(elmasks, 1) == nelements(mosaic)
    res = zeros(Int, set2index === nothing ? nsets(mosaic) : length(set2index), size(elmasks, 2))
    @inbounds for maskix in 1:size(elmasks, 2)
        nmasked_perset_mask!(view(res, :, maskix), mosaic, view(elmasks, :, maskix), set2index)
    end
    return res
end

# info for the set × mask overlap for the MaskedSetMosaic
struct MaskOverlap
    mask::Int       # index of the mask
    nmasked::Int    # number of masked elements
    nunmasked::Int  # number of unmasked elements
end

"""
`SetMosaic` with an elements mask (selection) on top.
Sets that are not overlapping with the mask are excluded(skipped) from `MaskedSetMosaic`.
Optionally the filtering can include testing for the minimal overlap significance P-value.

The tiles of non-overlapped sets are removed, the tiles that have identical membership
for all the masked sets are squashed into a single tile.
"""
mutable struct MaskedSetMosaic{T,S,M}
    original::SetMosaic{T,S}        # original mosaic
    elmasks::BitMatrix              # elements×masks
    elunionmask::BitVector          # all masked elements
    total_masked::Vector{Int}       # total number of masked elements for each mask

    # info for all non-empty overlaps with all masks
    # masks in MaskOverlap vector are ordered by mask indices
    set2masks::Dict{Int, Vector{MaskOverlap}}

    ix2mask::Vector{M}              # mask index to the ID of the mask
    mask2ix::Dict{M, Int}           # mask ID to index
end

function MaskedSetMosaic(mosaic::SetMosaic{T, S}, elmasks::AbstractMatrix{Bool},
                         set2masks::Dict{Int, Vector{MaskOverlap}},
                         mask_ids::Union{AbstractVector{M}, AbstractSet{M}, Nothing} = nothing
                        ) where {T,S,M}
    size(elmasks, 1) == nelements(mosaic) ||
        throw(ArgumentError("Elements mask length ($(size(elmasks, 1))) should match the number of elements ($(nelements(mosaic)))"))
    nmasks = size(elmasks, 2)
    (mask_ids === nothing) || (nmasks == length(mask_ids)) ||
        throw(ArgumentError("Number of masks ($nmasks) should match the number of mask IDs ($(length(masks_ids)))"))
    bit_elmasks = convert(BitMatrix, elmasks)
    bit_elunionmask = dropdims(any(bit_elmasks, dims=2), dims=2)
    @assert length(bit_elunionmask) == size(bit_elmasks, 1)
    MaskedSetMosaic(mosaic, bit_elmasks, bit_elunionmask,
                    dropdims(sum(elmasks, dims=1), dims=1), set2masks,
                    collect(mask_ids !== nothing ? mask_ids : 1:nmasks),
                    # FIXME use vector container if mask_ids === nothing or if ids are integers
                    Dict(id => ix for (ix, id) in enumerate(mask_ids !== nothing ? mask_ids : 1:nmasks)))
end

function mask(mosaic::SetMosaic, elmasks::AbstractMatrix{Bool};
              mask_ids::Union{AbstractVector, AbstractSet, Nothing} = nothing,
              min_nmasked::Integer=1, max_setsize::Union{Integer, Nothing} = nothing,
              max_overlap_logpvalue::Float64=0.0 # 0.0 would accept any overlap (as log(Fisher Exact Test P-value))
)
    size(elmasks, 1) == nelements(mosaic) ||
        throw(ArgumentError("Elements mask rows ($(size(elmasks, 1))) should match the number of elements ($(nelements(mosaic)))"))

    # calculate overlap sizes for each set and each mask
    nmasked_orgsets = nmasked_perset(mosaic, elmasks)

    # get the sets that overlap with the masked elements and with at least max_overlap_logpvalue significance
    set2masks = Dict{Int, Vector{MaskOverlap}}()
    emptyoverlap = Vector{MaskOverlap}()
    ntotal = nelements(mosaic)
    @inbounds for maskix in axes(nmasked_orgsets, 2)
        ntotal_masked = sum(view(elmasks, :, maskix))
        orgsets_mask = view(nmasked_orgsets, :, maskix)
        for (setix, nmasked) in enumerate(orgsets_mask)
            (nmasked < min_nmasked) && continue # skip small overlap
            @inbounds nset = setsize(mosaic, setix)
            (max_setsize !== nothing) && (nset > max_setsize) && continue # skip very big and generic sets
            overlap_pvalue = logpvalue(nmasked, nset, ntotal_masked, ntotal)
            (overlap_pvalue > max_overlap_logpvalue) && continue # skip non-signif overlaps
            setmasks = get!(() -> Vector{MaskOverlap}(), set2masks, setix)
            push!(setmasks, MaskOverlap(maskix, nmasked, setsize(mosaic, setix) - nmasked))
        end
    end

    return MaskedSetMosaic(mosaic, elmasks, set2masks, mask_ids)
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
         mask_ids = keys(elmasks), kwargs...)
end

unmask(mosaic::MaskedSetMosaic) = mosaic.original

nelements(mosaic::MaskedSetMosaic) = nelements(mosaic.original)
nsets(mosaic::MaskedSetMosaic) = length(mosaic.set2masks) # only sets overlapping with masks
nsets(mosaic::MaskedSetMosaic, maskix::Int) = # number of sets overlapping with given mask
    mapreduce(olaps -> any(olap -> olap.mask == maskix), sum, values(mosaic.set2masks), init=0)
nmasks(mosaic::MaskedSetMosaic) = length(mosaic.total_masked)
nmasked(mosaic::MaskedSetMosaic, maskix::Int) = mosaic.total_masked[maskix]
nunmasked(mosaic::MaskedSetMosaic, maskix::Int) = nelements(mosaic) - mosaic.total_masked[maskix]

# copy everything, except the original mosaic (leave the reference to the same object)
Base.copy(mosaic::MaskedSetMosaic) =
    MaskedSetMosaic(mosaic.original, copy(mosaic.elmasks), copy(mosaic.elunionmask),
                    copy(mosaic.total_masked),
                    deepcopy(mosaic.set2masks),
                    copy(mosaic.ix2mask), copy(mosaic.mask2ix))
