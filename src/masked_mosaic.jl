# the number of masked elements in each set of mosaic
function nmasked_perset(mosaic::SetMosaic, elmasks::AbstractVector)
    res = zeros(Int, nsets(mosaic), length(elmasks))
    @inbounds for (maskix, elmask) in enumerate(elmasks)
        @assert length(elmask) == nelements(mosaic)
        @inbounds for elix in eachindex(elmask)
            elmask[elix] || continue
            for setix in view(mosaic.setXelm, :, elix)
                res[setix, maskix] += 1
            end
        end
    end
    return res
end

# the number of masked elements in each set of mosaic
function nmasked_perset(mosaic::SetMosaic, elmasks::AbstractMatrix{Bool})
    @assert size(elmasks, 1) == nelements(mosaic)
    res = zeros(Int, nsets(mosaic), size(elmasks, 2))
    @inbounds for elix in 1:size(elmasks,1), maskix in 1:size(elmasks, 2)
        elmasks[elix, maskix] || continue
        for setix in view(mosaic.setXelm, :, elix)
            res[setix, maskix] += 1
        end
    end
    return res
end

"""
`SetMosaic` with an elements mask (selection) on top.
Sets that are not overlapping with the mask are excluded(skipped) from `MaskedSetMosaic`.
Optionally the filtering can include testing for the minimal overlap significance P-value.

The tiles of non-overlapped sets are removed, the tiles that have identical membership
for all the masked sets are squashed into a single tile.
"""
mutable struct MaskedSetMosaic{T,S}
    original::SetMosaic{T,S}    # original mosaic
    elmasks::BitMatrix          # elements masks
    setixs::Vector{Int}         # original indices of the included sets

    total_masked::Vector{Int}       # total number of masked elements
    nmasked_perset::Matrix{Int}     # number of masked elements in a set
    nunmasked_perset::Matrix{Int}   # number of unmasked elements in a set

    function MaskedSetMosaic(original::SetMosaic{T, S}, elmasks::BitMatrix,
                             setixs::Vector{Int}, total_masked::Vector{Int},
                             nmasked_perset::Matrix{Int}, nunmasked_perset::Matrix{Int}
    ) where {T,S}
        new{T,S}(original, elmasks, setixs,
                 total_masked, nmasked_perset, nunmasked_perset)
    end

    function MaskedSetMosaic(mosaic::SetMosaic{T,S}, elmasks::AbstractMatrix{Bool},
                             setixs::Vector{Int}, nmasked_perset::Matrix{Int}) where {T, S}
        size(elmasks, 1) == nelements(mosaic) ||
            throw(ArgumentError("Elements mask length ($(length(elmasks))) should match the number of elements ($(nelements(mosaic)))"))
        @assert nsets(mosaic) == size(nmasked_perset, 1)
        @assert size(elmasks, 2) == size(nmasked_perset, 2)
        # FIXME check setixs

        @inbounds nunmasked_newsets = [mosaic.set_sizes[setix] - nmasked_perset[setix, maskix] for setix in setixs, maskix in 1:size(nmasked_perset, 2)]
        new{T,S}(mosaic, elmasks, setixs, squeeze(sum(elmasks, 1), 1),
                 nmasked_perset[setixs, :], nunmasked_newsets)
    end

    MaskedSetMosaic(mosaic::SetMosaic{T, S},
                    elmasks::AbstractMatrix{Bool},
                    setixs::Vector{Int}) where {T,S} =
        MaskedSetMosaic(mosaic, elmasks, setixs, nmasked_perset(mosaic, elmasks))
end

# FIXME not optimal
setid2ix(mosaic::MaskedSetMosaic{T,S}, set::S) where {T,S} =
    searchsortedfirst(mosaic.setixs, mosaic.original.set2ix[set])

function mask(mosaic::SetMosaic, elmasks::AbstractMatrix{Bool};
              max_overlap_logpvalue::Float64 = 0.0 # 0.0 would accept any overlap (as log(Fisher Exact Test P-value))
)
    size(elmasks, 1) == nelements(mosaic) ||
        throw(ArgumentError("Elements mask rows ($(size(elmasks, 1))) should match the number of elements ($(nelements(mosaic)))"))

    # get the sets that overlap with the mask elements and with at least max_overlap_logpvalue significance
    nmasked_orgsets = nmasked_perset(mosaic, elmasks)
    ntotal = nelements(mosaic)
    nmasked = squeeze(sum(elmasks, 1), 1)
    org_setixs = sizehint!(Vector{Int}(), nsets(mosaic))
    for org_setix in 1:size(nmasked_orgsets, 1)
        for maskix in 1:size(nmasked_orgsets, 2)
            @inbounds nmasked_orgset = nmasked_orgsets[org_setix, maskix]
            (nmasked_orgset == 0) && continue
            @inbounds overlap_pvalue = logpvalue(nmasked_orgset, setsize(mosaic, org_setix), nmasked[maskix], ntotal)
            if overlap_pvalue <= max_overlap_logpvalue
                # add the set, stop the loop
                push!(org_setixs, org_setix)
                break
            end
        end
    end

    return MaskedSetMosaic(mosaic, elmasks, org_setixs, nmasked_orgsets)
end

mask(mosaic::SetMosaic{T,S}, elmasks::Vector{Set{T}}; max_overlap_logpvalue::Real = 0.0) where {T,S} =
    mask(mosaic, Bool[in(e, elmask) for e in mosaic.ix2elm, elmask in elmasks],
         max_overlap_logpvalue=max_overlap_logpvalue)
unmask(mosaic::MaskedSetMosaic) = mosaic.original

function setsize(mosaic::MaskedSetMosaic, setix::Integer)
    @inbounds org_setix = mosaic.setixs[setix]
    return setsize(mosaic.original, org_setix)
end

nelements(mosaic::MaskedSetMosaic) = nelements(mosaic.original)
nsets(mosaic::MaskedSetMosaic) = length(mosaic.setixs)
nmasks(mosaic::MaskedSetMosaic) = length(mosaic.total_masked)
nmasked(mosaic::MaskedSetMosaic, maskix::Int) = mosaic.total_masked[maskix]
nunmasked(mosaic::MaskedSetMosaic, maskix::Int) = nelements(mosaic) - mosaic.total_masked[maskix]
nmasked_perset(mosaic::MaskedSetMosaic, maskix::Int) = view(mosaic.nmasked_perset, :, maskix)
nunmasked_perset(mosaic::MaskedSetMosaic, maskix::Int) = view(mosaic.nunmasked_perset, :, maskix)
nmasked(mosaic::MaskedSetMosaic{T,S}, set::S, maskix::Int) where {T,S} =
    (setix = setid2ix(mosaic, set); setix > 0 ? mosaic.nmasked_perset[setix, maskix] : 0)
nunmasked(mosaic::MaskedSetMosaic{T,S}, set::S, maskix::Int) where {T,S} =
    (setix = setid2ix(mosaic, set); setix > 0 ? mosaic.nunmasked_perset[setix, maskix] : setsize(mosaic.original, mosaic.original.set2ix[set]))

function Base.copy(mosaic::MaskedSetMosaic)
    # copy everything, except the original mosaic (leave the reference to the same object)
    MaskedSetMosaic(mosaic.original, copy(mosaic.elmasks), copy(mosaic.setixs),
                    mosaic.total_masked,
                    copy(mosaic.nmasked_perset), copy(mosaic.nunmasked_perset))
end

"""
Exclude `setmask` sets from the `mosaic` and update the set of its active tiles.
"""
function Base.filter!(mosaic::MaskedSetMosaic, setmask::Union{Vector{Bool},BitVector})
    nsets(mosaic) == length(setmask) ||
        throw(ArgumentError("Mask length ($(length(setmask))) does not match the number of sets in mosaic ($(nsets(mosaic)))"))
    mosaic.setixs = mosaic.setixs[setmask]
    mosaic.nmasked_perset = mosaic.nmasked_perset[setmask, :]
    mosaic.nunmasked_perset = mosaic.nunmasked_perset[setmask, :]
    return mosaic
end
