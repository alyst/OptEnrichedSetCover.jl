# the number of masked elements in each set of mosaic
function nmasked_perset_mask!(nmasked::AbstractVector{Int}, mosaic::SetMosaic, elmask::AbstractVector{Bool},
                              set2index::Union{Void,Dict{Int,Int}} = nothing)
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
                        set2index::Union{Void,Dict{Int,Int}} = nothing)
    res = zeros(Int, set2index === nothing ? nsets(mosaic) : length(set2index), length(elmasks))
    @inbounds for (maskix, elmask) in enumerate(elmasks)
        @assert length(elmask) == nelements(mosaic)
        nmasked_perset_mask!(view(res, :, maskix), mosaic, elmask, set2index)
    end
    return res
end

# the number of masked elements in each set of mosaic
function nmasked_perset(mosaic::SetMosaic, elmasks::AbstractMatrix{Bool},
                        set2index::Union{Void,Dict{Int,Int}} = nothing)
    @assert size(elmasks, 1) == nelements(mosaic)
    res = zeros(Int, set2index === nothing ? nsets(mosaic) : length(set2index), size(elmasks, 2))
    @inbounds for maskix in 1:size(elmasks, 2)
        nmasked_perset_mask!(view(res, :, maskix), mosaic, view(elmasks, :, maskix), set2index)
    end
    return res
end

# info for the set from the MaskedSetMosaic
struct MaskedSet
    mask::Int       # index of the mask
    set::Int        # index of the set in the SetMosaic
    nmasked::Int    # number of masked elements
    nunmasked::Int  # number of unmasked elements
end

setsize(s::MaskedSet) = s.nmasked + s.nunmasked

"""
`SetMosaic` with an elements mask (selection) on top.
Sets that are not overlapping with the mask are excluded(skipped) from `MaskedSetMosaic`.
Optionally the filtering can include testing for the minimal overlap significance P-value.

The tiles of non-overlapped sets are removed, the tiles that have identical membership
for all the masked sets are squashed into a single tile.
"""
mutable struct MaskedSetMosaic{T,S}
    original::SetMosaic{T,S}        # original mosaic
    elmasks::BitMatrix              # elements masks
    total_masked::Vector{Int}       # total number of masked elements for each mask
    maskedsets::Vector{MaskedSet}   # per set mask info (only the sets intersecting with each mask included)
    orig2masked::Dict{Int, Set{Int}}# map from the original index to the indices of masked sets
end

function MaskedSetMosaic(mosaic::SetMosaic{T, S}, elmasks::AbstractMatrix{Bool},
                         maskedsets::Vector{MaskedSet}) where {T,S}
    size(elmasks, 1) == nelements(mosaic) ||
        throw(ArgumentError("Elements mask length ($(length(elmasks))) should match the number of elements ($(nelements(mosaic)))"))
    nmasks = size(elmasks, 2)
    orig2masked = Dict{Int,Set{Int}}()
    for (i, mset) in enumerate(maskedsets)
        1 <= mset.mask <= nmasks || throw(ArgumentError("Incorrect mask index ($(mset.mask)) for set #$i"))
        1 <= mset.set <= nsets(mosaic) || throw(ArgumentError("Incorrect set index ($(mset.set) for set #$i"))
        maskixs = get!(() -> Set{Int}(), orig2masked, mset.set)
        push!(maskixs, i)
    end
    MaskedSetMosaic(mosaic, convert(BitMatrix, elmasks),
                    squeeze(sum(elmasks, 1), 1), maskedsets, orig2masked)
end

# FIXME not optimal
setid2ix(mosaic::MaskedSetMosaic{T,S}, set::S) where {T,S} =
    searchsortedfirst(mosaic.setixs, mosaic.original.set2ix[set])

function mask(mosaic::SetMosaic, elmasks::AbstractMatrix{Bool};
              min_nmasked::Integer=2,
              max_overlap_logpvalue::Float64=0.0 # 0.0 would accept any overlap (as log(Fisher Exact Test P-value))
)
    size(elmasks, 1) == nelements(mosaic) ||
        throw(ArgumentError("Elements mask rows ($(size(elmasks, 1))) should match the number of elements ($(nelements(mosaic)))"))

    # get the sets that overlap with the mask elements and with at least max_overlap_logpvalue significance
    nmasked_orgsets = nmasked_perset(mosaic, elmasks)
    ntotal = nelements(mosaic)
    maskedsets = Vector{MaskedSet}()

    @inbounds for maskix in 1:size(nmasked_orgsets, 2)
        ntotal_masked = sum(view(elmasks, :, maskix))
        orgsets_mask = view(nmasked_orgsets, :, maskix)
        for (setix, nmasked) in enumerate(orgsets_mask)
            (nmasked < min_nmasked) && continue
            nset = setsize(mosaic, setix)
            @inbounds overlap_pvalue = logpvalue(nmasked, nset, ntotal_masked, ntotal)
            if overlap_pvalue <= max_overlap_logpvalue
                push!(maskedsets, MaskedSet(maskix, setix, nmasked, nset - nmasked))
            end
        end
    end

    return MaskedSetMosaic(mosaic, elmasks, maskedsets)
end

mask(mosaic::SetMosaic{T,S}, elmasks::Vector{Set{T}};
     min_nmasked::Integer=2, max_overlap_logpvalue::Real=0.0) where {T,S} =
    mask(mosaic, Bool[in(e, elmask) for e in mosaic.ix2elm, elmask in elmasks],
         min_nmasked=min_nmasked, max_overlap_logpvalue=max_overlap_logpvalue)
unmask(mosaic::MaskedSetMosaic) = mosaic.original

nelements(mosaic::MaskedSetMosaic) = nelements(mosaic.original)
nsets(mosaic::MaskedSetMosaic) = length(mosaic.maskedsets)
nsets(mosaic::MaskedSetMosaic, maskix::Int) = sum(ms -> ms.mask == maskix, mosaic.maskedsets)
nmasks(mosaic::MaskedSetMosaic) = length(mosaic.total_masked)
nmasked(mosaic::MaskedSetMosaic, maskix::Int) = mosaic.total_masked[maskix]
nunmasked(mosaic::MaskedSetMosaic, maskix::Int) = nelements(mosaic) - mosaic.total_masked[maskix]
maskedset(mosaic::MaskedSetMosaic, ix::Int) = mosaic.maskedsets[ix]
nmasked(mosaic::MaskedSetMosaic{T,S}, set::S, maskix::Int) where {T,S} =
    (setix = setid2ix(mosaic, set); setix > 0 ? mosaic.nmasked_perset[setix, maskix] : 0)
nunmasked(mosaic::MaskedSetMosaic{T,S}, set::S, maskix::Int) where {T,S} =
    (setix = setid2ix(mosaic, set); setix > 0 ? mosaic.nunmasked_perset[setix, maskix] : setsize(mosaic.original, mosaic.original.set2ix[set]))

function Base.copy(mosaic::MaskedSetMosaic)
    # copy everything, except the original mosaic (leave the reference to the same object)
    MaskedSetMosaic(mosaic.original, copy(mosaic.elmasks), copy(mosaic.total_masked),
                    copy(mosaic.maskedsets), copy(mosaic.orig2masked))
end

"""
Exclude `setmask` sets from the `mosaic` and update the set of its active tiles.
"""
function Base.filter!(mosaic::MaskedSetMosaic, setmask::AbstractVector{Bool})
    nsets(mosaic) == length(setmask) ||
        throw(ArgumentError("Mask length ($(length(setmask))) does not match the number of sets in mosaic ($(nsets(mosaic)))"))
    mosaic.maskedsets = mosaic.maskedsets[setmask]
    return mosaic
end
