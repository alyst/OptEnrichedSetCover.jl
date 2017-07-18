function _encode_elements{T}(all_elms::Set{T})
    ix2elm = collect(all_elms)
    sort!(ix2elm)
    # map from element to its index
    elm2ix = Dict{T,Int}()
    sizehint!(elm2ix, length(ix2elm))
    for i in 1:length(ix2elm)
        elm2ix[ix2elm[i]] = i
    end
    ix2elm, elm2ix
end

function _prepare_tiles{T}(sets, elm2ix::Dict{T, Int})
    if isempty(sets)
        # no sets -- no tiles
        return SparseMaskMatrix(0, length(elm2ix)), SparseMaskMatrix(length(elm2ix), 0),
               SparseMaskMatrix(length(elm2ix), 0), SparseMaskMatrix(0, 0)
    end

    # build set-X-element membership matrix
    elmXset = fill(false, length(elm2ix), length(sets))
    for (i, s) in enumerate(sets)
        set_elms = view(elmXset, :, i)
        for e in s
            @inbounds set_elms[elm2ix[e]] = true
        end
    end
    setXelm = elmXset'
    # squash element-X-set matrix (duplicate rows) into membership-mask-to-tile dict
    sets2elms = Dict{Vector{Int}, Vector{Int}}()
    sizehint!(sets2elms, length(elm2ix))
    for i in 1:size(setXelm, 2)
        @inbounds set_ixs = find(view(setXelm, :, i))
        push!(get!(() -> Vector{Int}(), sets2elms, set_ixs), i)
    end

    # build tile-X-set and element-X-tile membership matrices
    set2tile_ixs = [Vector{Int}() for _ in 1:length(sets)]
    tile_elm_ranges = Vector{Int}()
    elm_ixs = Vector{Int}()
    for (set_ixs, tile_elm_ixs) in sets2elms
        push!(tile_elm_ranges, length(elm_ixs)+1)
        tile_ix = length(tile_elm_ranges)
        append!(elm_ixs, tile_elm_ixs)

        for set_ix in set_ixs
            push!(set2tile_ixs[set_ix], tile_ix)
        end
    end
    push!(tile_elm_ranges, length(elm_ixs)+1)

    set_tile_ranges = Vector{Int}()
    tile_ixs = Vector{Int}()
    for (set_ix, set_tile_ixs) in enumerate(set2tile_ixs)
        push!(set_tile_ranges, length(tile_ixs)+1)
        append!(tile_ixs, set_tile_ixs)
    end
    push!(set_tile_ranges, length(tile_ixs)+1)

    SparseMaskMatrix(setXelm),
    SparseMaskMatrix(elmXset),
    SparseMaskMatrix(length(elm2ix), length(tile_elm_ranges)-1, tile_elm_ranges, elm_ixs), # elmXtile
    SparseMaskMatrix(length(tile_elm_ranges)-1, length(sets), set_tile_ranges, tile_ixs) # tileXset
end

# size of the sets intersection
# the sets are represented by the sorted sets of their tiles
function _isect_size(set1_tiles::StridedVector{Int},
                     set2_tiles::StridedVector{Int},
                     tile_sizes::Vector{Int})
    (length(set1_tiles)==0 || length(set2_tiles)==0) && return 0
    res = 0
    i = 1
    tile1 = set1_tiles[i] # current tile of the 1st set
    j = 1
    tile2 = set2_tiles[j] # current tile of the 2nd set
    @inbounds while true
        if tile1 < tile2
            i += 1
            (i <= length(set1_tiles)) || break
            tile1 = set1_tiles[i]
        elseif tile1 > tile2
            j += 1
            (j <= length(set2_tiles)) || break
            tile2 = set2_tiles[j]
        else # intersecting tile, adjust the overlap
            res += tile_sizes[tile1]
            i += 1
            j += 1
            (i <= length(set1_tiles) && j <= length(set2_tiles)) || break
            tile1 = set1_tiles[i]
            tile2 = set2_tiles[j]
        end
    end
    return res
end

function _setXset_scores(tileXset::SparseMaskMatrix, total_size::Int, set_sizes::Vector{Int}, tile_sizes::Vector{Int})
    nsets = size(tileXset, 2)
    res = Matrix{Float64}(nsets, nsets)
    @inbounds for set1_ix in 1:nsets
        res[set1_ix, set1_ix] = 0.0
        (set1_ix == nsets) && break
        set1_tiles = view(tileXset, :, set1_ix)
        for set2_ix in (set1_ix+1):nsets
            set2_tiles = view(tileXset, :, set2_ix)
            # one-sided Fisher's P-value, right tail
            res[set2_ix, set1_ix] =
                logpvalue(set_sizes[set1_ix], set_sizes[set2_ix], total_size,
                          _isect_size(set1_tiles, set2_tiles, tile_sizes))
        end
    end
    # symmetrize
    @inbounds for set2_ix in 1:nsets
        for set1_ix in (set2_ix+1):nsets
            res[set2_ix, set1_ix] = res[set1_ix, set2_ix]
        end
    end
    return res
end

# unlike Base.union() does not use recursion
function _union{T}(::Type{T}, sets)
    res = Set{T}()
    for set in sets
        union!(res, set)
    end
    return res
end

"""
A collection of (potentially overlapping) sets as
a "mosaic" of non-overlapping "tiles".

* `T` type of elements
* `S` type of set keys
"""
@compat struct SetMosaic{T,S}
    ix2elm::Vector{T}       # element index to element
    elm2ix::Dict{T, Int}    # element to its index

    ix2set::Vector{S}       # set index to ID
    set2ix::Dict{S, Int}    # set ID to index

    set_sizes::Vector{Int}

    setXelm::SparseMaskMatrix  # set×element membership
    elmXset::SparseMaskMatrix  # element×set membership (transpose of setXelm)
    elmXtile::SparseMaskMatrix  # element-to-tile mask
    tileXset::SparseMaskMatrix  # rows = tiles, cols = sets from collection, true if tile is a subset of a set

    setXset_scores::Matrix{Float64}

    """
    Construct `SetMosaic` for a given nameless sets collection.
    """
    function (::Type{SetMosaic}){T}(sets::Vector{Set{T}}, all_elms::Set{T} = _union(T, sets))
        ix2elm, elm2ix = _encode_elements(all_elms)
        setXelm, elmXset, elmXtile, tileXset = _prepare_tiles(sets, elm2ix)
        tile_sizes = Int[length(view(elmXtile, :, tile_ix)) for tile_ix in 1:size(elmXtile, 2)]
        set_sizes = Int[length(set) for set in sets]
        new{T, Int}(ix2elm, elm2ix,
                    collect(eachindex(sets)), Dict([Pair(i,i) for i in eachindex(sets)]),
                    set_sizes,
                    setXelm, elmXset, elmXtile, tileXset,
                    _setXset_scores(tileXset, length(all_elms), set_sizes, tile_sizes))
    end

    """
    Constructs `SetMosaic` for a given named sets collection.
    """
    function (::Type{SetMosaic}){T,S}(sets::Dict{S, Set{T}}, all_elms::Set{T} = _union(T, values(sets)))
        ix2elm, elm2ix = _encode_elements(all_elms)
        setXelm, elmXset, elmXtile, tileXset = _prepare_tiles(values(sets), elm2ix)
        tile_sizes = Int[length(view(elmXtile, :, tile_ix)) for tile_ix in 1:size(elmXtile, 2)]
        set_sizes = Int[length(set) for set in values(sets)]
        new{T, S}(ix2elm, elm2ix,
                  collect(keys(sets)), Dict([Pair(s, i) for (i, s) in enumerate(keys(sets))]),
                  set_sizes,
                  setXelm, elmXset, elmXtile, tileXset,
                  _setXset_scores(tileXset, length(all_elms), set_sizes, tile_sizes))
    end
end

nelements(mosaic::SetMosaic) = length(mosaic.ix2elm)
ntiles(mosaic::SetMosaic) = size(mosaic.tileXset, 1)
nsets(mosaic::SetMosaic) = size(mosaic.tileXset, 2)

tile(mosaic::SetMosaic, tile_ix::Integer) = view(mosaic.elmXtile, :, tile_ix)
set(mosaic::SetMosaic, set_ix::Integer) = view(mosaic.elmXset, :, set_ix)

setsize(mosaic::SetMosaic, set_ix::Integer) = mosaic.set_sizes[set_ix]

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
type MaskedSetMosaic{T,S}
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

    function MaskedSetMosaic(mosaic::SetMosaic{T,S}, elmasks::Union{BitMatrix, Matrix{Bool}},
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

    function MaskedSetMosaic(mosaic::SetMosaic{T,S}, elmasks::Union{BitMatrix, Matrix{Bool}},
                             max_overlap_logpvalue::Float64 = 0.0 # 0.0 would accept any overlap (as log(Fisher Exact Test P-value))
    ) where {T,S}
        size(elmasks, 1) == nelements(mosaic) ||
            throw(ArgumentError("Elements mask rows ($(size(elmasks, 1))) should match the number of elements ($(nelements(mosaic)))"))
        (max_overlap_logpvalue <= 0.0) || throw(ArgumentError("Maximal overlap log(P-value) must be ≤0, found $max_overlap_logpvalue"))

        # get the sets that overlap with the mask elements and with at least max_overlap_logpvalue significance
        nmasked_orgsets = nmasked_perset(mosaic, elmasks)
        ntotal = nelements(mosaic)
        nmasked = squeeze(sum(elmasks, 1), 1)
        org_setixs = sizehint!(Vector{Int}(), nsets(mosaic))
        for org_setix in 1:size(nmasked_orgsets, 1)
            for maskix in 1:size(nmasked_orgsets, 2)
                @inbounds nmasked_orgset = nmasked_orgsets[org_setix, maskix]
                (nmasked_orgset == 0) && continue
                @inbounds overlap_pvalue = logpvalue(setsize(mosaic, org_setix), nmasked[maskix], ntotal, nmasked_orgset)
                if overlap_pvalue <= max_overlap_logpvalue
                    # add the set, stop the loop
                    push!(org_setixs, org_setix)
                    break
                end
            end
        end

        return MaskedSetMosaic(mosaic, elmasks, org_setixs, nmasked_orgsets)
    end

    MaskedSetMosaic(mosaic::SetMosaic{T, S},
                    elmasks::Union{BitMatrix, Matrix{Bool}},
                    setixs::Vector{Int}) where {T,S} =
        MaskedSetMosaic(mosaic, elmasks, setixs, nmasked_perset(mosaic, elmasks))
end

# FIXME not optimal
setid2ix{T,S}(mosaic::MaskedSetMosaic{T,S}, set::S) = searchsortedfirst(mosaic.setixs, mosaic.original.set2ix[set])

mask(mosaic::SetMosaic, mask::Union{BitMatrix,Matrix{Bool}}; max_overlap_logpvalue::Real = 0.0) =
    MaskedSetMosaic(mosaic, mask, max_overlap_logpvalue)
mask{T,S}(mosaic::SetMosaic{T,S}, sels::Vector{Set{T}}; max_overlap_logpvalue::Real = 0.0) =
    mask(mosaic, Bool[in(e, sel) for e in mosaic.ix2elm, sel in sels],
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
nmasked{T,S}(mosaic::MaskedSetMosaic{T,S}, set::S, maskix::Int) = (setix = setid2ix(mosaic, set); setix > 0 ? mosaic.nmasked_perset[setix, maskix] : 0)
nunmasked{T,S}(mosaic::MaskedSetMosaic{T,S}, set::S, maskix::Int) = (setix = setid2ix(mosaic, set); setix > 0 ? mosaic.nunmasked_perset[setix, maskix] : setsize(mosaic.original, mosaic.original.set2ix[set]))

function Base.copy{T,S}(mosaic::MaskedSetMosaic{T,S})
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
