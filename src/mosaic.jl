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
immutable SetMosaic{T,S}
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

"""
`SetMosaic` with an elements mask (selection) on top.
Sets that are not overlapping with the mask are excluded(skipped) from `MaskedSetMosaic`.
Optionally the filtering can include testing for the minimal overlap significance P-value.

The tiles of non-overlapped sets are removed, the tiles that have identical membership
for all the masked sets are squashed into a single tile.
"""
type MaskedSetMosaic{T,S}
    original::SetMosaic{T,S}    # original mosaic
    elmask::BitVector           # elements mask
    setixs::Vector{Int}         # original indices of the included sets

    tileXset::SparseMaskMatrix
    total_masked::Int               # total number of masked elements
    nmasked_pertile::Vector{Int}    # number of masked elements in a tile
    nunmasked_pertile::Vector{Int}  # number of unmasked elements in a tile
    nmasked_perset::Vector{Int}     # number of masked elements in a set
    nunmasked_perset::Vector{Int}   # number of unmasked elements in a set

    function MaskedSetMosaic(original::SetMosaic{T, S}, elmask::BitVector,
                             setixs::Vector{Int},
                             tileXset::SparseMaskMatrix,
                             total_masked::Number,
                             nmasked_pertile::Vector{Int}, nunmasked_pertile::Vector{Int},
                             nmasked_perset::Vector{Int}, nunmasked_perset::Vector{Int}
    )
        new(original, elmask, setixs, tileXset, total_masked,
            nmasked_pertile, nunmasked_pertile,
            nmasked_perset, nunmasked_perset)
    end

    function (::Type{MaskedSetMosaic}){T,S}(mosaic::SetMosaic{T, S}, elmask::Union{BitVector, Vector{Bool}},
                                            max_overlap_logpvalue::Float64 = 0.0 # 0.0 would accept any overlap (as log(Fisher Exact Test P-value))
    )
        length(elmask) == nelements(mosaic) ||
            throw(ArgumentError("Elements mask length ($(length(elmask))) should match the number of elements ($(nelements(mosaic)))"))
        (max_overlap_logpvalue <= 0.0) || throw(ArgumentError("Maximal overlap log(P-value) must be ≤0, found $max_overlap_logpvalue"))

        # get the sets that overlap with the mask elements and with at least max_overlap_logpvalue significance
        nmasked_orgsets = zeros(Int, nsets(mosaic))
        @inbounds for elix in eachindex(elmask)
            elmask[elix] || continue
            for org_setix in view(mosaic.setXelm, :, elix)
                nmasked_orgsets[org_setix] += 1
            end
        end
        ntotal = nelements(mosaic)
        nmasked = sum(elmask)
        org_setixs = sizehint!(Vector{Int}(), nsets(mosaic))
        for (org_setix, nmasked_orgset) in enumerate(nmasked_orgsets)
            @inbounds overlap_pvalue = logpvalue(setsize(mosaic, org_setix), nmasked, ntotal, nmasked_orgset)
            if (nmasked_orgset > 0) && (overlap_pvalue <= max_overlap_logpvalue)
                push!(org_setixs, org_setix)
            end
        end

        # detect and squash tiles that have the same membership pattern wrt the masked sets
        newsetXorgtile = fill(false, length(org_setixs), ntiles(mosaic))
        @inbounds for (new_setix, org_setix) in enumerate(org_setixs)
            newsetXorgtile[new_setix, view(mosaic.tileXset, :, org_setix)] = true
        end
        newsets2orgtiles = Dict{Vector{Int}, Vector{Int}}()
        sizehint!(newsets2orgtiles, size(newsetXorgtile, 2))
        for org_tileix in 1:size(newsetXorgtile, 2)
            @inbounds new_setixs = find(view(newsetXorgtile, :, org_tileix))
            isempty(new_setixs) && continue
            push!(get!(() -> Vector{Int}(), newsets2orgtiles, new_setixs), org_tileix)
        end

        # build newtile-to-newset membership matrix
        nmasked_newtiles = zeros(Int, length(newsets2orgtiles))
        nunmasked_newtiles = zeros(Int, length(newsets2orgtiles))
        newset2newtiles = [Vector{Int}() for _ in eachindex(org_setixs)]
        @inbounds for (new_tileix, (new_setixs, org_tileixs)) in enumerate(newsets2orgtiles)
            for new_setix in new_setixs
                push!(newset2newtiles[new_setix], new_tileix)
            end
            for org_tileix in org_tileixs
                orgtile_elms = view(mosaic.elmXtile, :, org_tileix)
                nmasked_orgtile = count(eix -> elmask[eix], orgtile_elms)
                nmasked_newtiles[new_tileix] += nmasked_orgtile
                nunmasked_newtiles[new_tileix] += length(orgtile_elms) - nmasked_orgtile
            end
        end
        newtileXnewset = SparseMaskMatrix(length(nmasked_newtiles), newset2newtiles)

        # calculate masked/unmasked elements for each set in the masked mosaic
        @inbounds nunmasked_newsets = [mosaic.set_sizes[org_setix] - nmasked_orgsets[org_setix] for org_setix in org_setixs]
        new{T,S}(mosaic, elmask, org_setixs, newtileXnewset,
                 sum(elmask), nmasked_newtiles, nunmasked_newtiles,
                 nmasked_orgsets[org_setixs], nunmasked_newsets)
    end
end

mask(mosaic::SetMosaic, mask::Union{BitVector,Vector{Bool}}; max_overlap_logpvalue::Real = 0.0) =
    MaskedSetMosaic(mosaic, mask, max_overlap_logpvalue)
mask{T,S}(mosaic::SetMosaic{T,S}, sel::Set{T}; max_overlap_logpvalue::Real = 0.0) =
    mask(mosaic, Bool[in(e, sel) for e in mosaic.ix2elm], max_overlap_logpvalue=max_overlap_logpvalue)
unmask(mosaic::MaskedSetMosaic) = mosaic.original

nelements(mosaic::MaskedSetMosaic) = nelements(mosaic.original)
ntiles(mosaic::MaskedSetMosaic) = size(mosaic.tileXset, 1)
nsets(mosaic::MaskedSetMosaic) = size(mosaic.tileXset, 2)
nmasked(mosaic::MaskedSetMosaic) = mosaic.total_masked
nunmasked(mosaic::MaskedSetMosaic) = nelements(mosaic) - mosaic.total_masked
nmasked_pertile(mosaic::MaskedSetMosaic) = mosaic.nmasked_pertile
nunmasked_pertile(mosaic::MaskedSetMosaic) = mosaic.nunmasked_pertile
nmasked_perset(mosaic::MaskedSetMosaic) = mosaic.nmasked_perset
nunmasked_perset(mosaic::MaskedSetMosaic) = mosaic.nunmasked_perset

function Base.copy{T,S}(mosaic::MaskedSetMosaic{T,S})
    # copy everything, except the original mosaic (leave the reference to the same object)
    MaskedSetMosaic{T,S}(mosaic.original, copy(mosaic.elmask), copy(mosaic.setixs),
                    copy(mosaic.tileXset), mosaic.total_masked,
                    copy(mosaic.nmasked_pertile), copy(mosaic.nunmasked_pertile),
                    copy(mosaic.nmasked_perset), copy(mosaic.nunmasked_perset))
end

"""
Exclude `setmask` sets from the `mosaic` and update the set of its active tiles.
"""
function Base.filter!(mosaic::MaskedSetMosaic, setmask::Union{Vector{Bool},BitVector})
    nsets(mosaic) == length(setmask) ||
        throw(ArgumentError("Mask length ($(length(setmask))) does not match the number of sets in mosaic ($(nsets(mosaic)))"))
    mosaic.setixs = mosaic.setixs[setmask]
    tileXset = mosaic.tileXset[:, setmask]
    tile_mask = fill(false, ntiles(mosaic))
    tile_mask[tileXset.rowval] = true
    old2new_tileixs = cumsum(tile_mask)
    old_tileixs = find(tile_mask)
    mosaic.tileXset = SparseMaskMatrix(length(old_tileixs), tileXset.n,
                                       tileXset.colptr, old2new_tileixs[tileXset.rowval]) # remove unused tiles
    mosaic.nmasked_pertile = mosaic.nmasked_pertile[old_tileixs]
    mosaic.nunmasked_pertile = mosaic.nunmasked_pertile[old_tileixs]
    mosaic.nmasked_perset = mosaic.nmasked_perset[setmask]
    mosaic.nunmasked_perset = mosaic.nunmasked_perset[setmask]
    return mosaic
end
