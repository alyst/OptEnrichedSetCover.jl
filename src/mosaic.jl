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
        return SparseMaskMatrix(0, length(elm2ix)), SparseMaskMatrix(length(elm2ix), 0), SparseMaskMatrix(0, 0)
    end

    # build set-X-element membership matrix
    setXelm = fill(false, length(sets), length(elm2ix))
    for (i, s) in enumerate(sets)
        for e in s
            setXelm[i, elm2ix[e]] = true
        end
    end
    # squash element-X-set matrix (duplicate rows) into membership-mask-to-tile dict
    membership2tile = Dict{Vector{Int}, Vector{Int}}()
    sizehint!(membership2tile, length(elm2ix))
    for i in 1:size(setXelm, 2)
        set_ixs = find(setXelm[:, i])
        if haskey(membership2tile, set_ixs)
            push!(membership2tile[set_ixs], i)
        else
            membership2tile[set_ixs] = fill(i, 1)
        end
    end

    # build tile-X-set and element-X-tile membership matrices
    set2tile_ixs = [Vector{Int}() for _ in 1:length(sets)]
    tile_elm_ranges = Vector{Int}()
    elm_ixs = Vector{Int}()
    for (set_ixs, tile_elm_ixs) in membership2tile
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
    SparseMaskMatrix(length(elm2ix), length(tile_elm_ranges)-1, tile_elm_ranges, elm_ixs), # elmXtile
    SparseMaskMatrix(length(tile_elm_ranges)-1, length(sets), set_tile_ranges, tile_ixs) # tileXset
end

# size of the sets intersection
# the sets are represented by the sorted sets of their tiles
function _isect_size(set1_ixs::StridedVector{Int}, set2_ixs::StridedVector{Int}, tile_sizes::Vector{Int})
    j = 1
    i = 1
    res = 0
    set1_ix = 0
    set2_ix = 0
    @inbounds while i <= length(set1_ixs) && j <= length(set2_ixs)
        if set1_ix==0
            set1_ix = set1_ixs[i]
        end
        if set2_ix==0
            set2_ix = set2_ixs[j]
        end
        if set1_ix < set2_ix
            i += 1
            set1_ix = 0
        elseif set1_ix > set2_ix
            j += 1
            set2_ix = 0
        else
            res += tile_sizes[set1_ix]
            i += 1
            j += 1
            set1_ix = set2_ix = 0
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
            # one-sided Fisher's P-value
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

    setXelm::SparseMaskMatrix  # set-to-element mask
    elmXtile::SparseMaskMatrix  # element-to-tile mask
    tileXset::SparseMaskMatrix  # rows = tiles, cols = sets from collection, true if tile is a subset of a set

    setXset_scores::Matrix{Float64}

    """
    Construct `SetMosaic` for a given nameless sets collection.
    """
    function (::Type{SetMosaic}){T}(sets::Vector{Set{T}}, all_elms::Set{T} = _union(T, sets))
        ix2elm, elm2ix = _encode_elements(all_elms)
        elmXset, elmXtile, tileXset = _prepare_tiles(sets, elm2ix)
        tile_sizes = Int[length(view(elmXtile, :, tile_ix)) for tile_ix in 1:size(elmXtile, 2)]
        set_sizes = Int[length(set) for set in sets]
        new{T, Int}(ix2elm, elm2ix,
                    collect(eachindex(sets)), Dict(Pair{Int,Int}[Pair(i,i) for i in eachindex(sets)]),
                    set_sizes,
                    elmXset, elmXtile, tileXset, _setXset_scores(tileXset, length(all_elms), set_sizes, tile_sizes))
    end

    """
    Constructs `SetMosaic` for a given named sets collection.
    """
    function (::Type{SetMosaic}){T,S}(sets::Dict{S, Set{T}}, all_elms::Set{T} = _union(T, values(sets)))
        ix2elm, elm2ix = _encode_elements(all_elms)
        elmXset, elmXtile, tileXset = _prepare_tiles(values(sets), elm2ix)
        tile_sizes = Int[length(view(elmXtile, :, tile_ix)) for tile_ix in 1:size(elmXtile, 2)]
        set_sizes = Int[length(set) for set in values(sets)]
        new{T, S}(ix2elm, elm2ix,
                  collect(keys(sets)), Dict{S, Int}(Pair{S,Int}[Pair(s, i) for (i, s) in enumerate(keys(sets))]),
                  set_sizes,
                  elmXset, elmXtile, tileXset, _setXset_scores(tileXset, length(all_elms), set_sizes, tile_sizes))
    end
end

nelements(mosaic::SetMosaic) = length(mosaic.ix2elm)
ntiles(mosaic::SetMosaic) = size(mosaic.tileXset, 1)
nsets(mosaic::SetMosaic) = size(mosaic.tileXset, 2)

tile(mosaic::SetMosaic, tile_ix::Integer) = view(mosaic.elmXtile, :, tile_ix)
set(mosaic::SetMosaic, set_ix::Integer) = view(mosaic.elmXset, :, set_ix)

"""
`SetMosaic` with an elements mask (selection) on top.
Sets that are not overlapping with the mask are excluded(skipped) from `MaskedSetMosaic`.

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

    function (::Type{MaskedSetMosaic}){T,S}(mosaic::SetMosaic{T, S}, elmask::Union{BitVector, Vector{Bool}})
        length(elmask) == nelements(mosaic) ||
            throw(ArgumentError("Elements mask length ($(length(elmask))) should match the number of elements ($(nelements(mosaic)))"))
        #println("elmask=", elmask)
        # get the sets that overlap with the elements mask
        setmask = fill(false, nsets(mosaic))
        @inbounds for i in eachindex(elmask)
            elmask[i] && (setmask[view(mosaic.setXelm, :, i)] = true)
        end
        orig_setixs = find(setmask)
        #println("setmask=$setmask")
        # detect and squash tiles that have the same membership pattern wrt the masked sets
        subsetXtile = fill(false, length(orig_setixs), ntiles(mosaic))
        for (new_setix, orig_setix) in enumerate(orig_setixs)
            subsetXtile[new_setix, view(mosaic.tileXset, :, orig_setix)] = true
        end
        #println("subsetXtile=$subsetXtile")
        membership2tile = Dict{Vector{Int}, Vector{Int}}()
        sizehint!(membership2tile, size(subsetXtile, 2))
        @inbounds for i in 1:size(subsetXtile, 2)
            set_ixs = find(subsetXtile[:, i])
            if haskey(membership2tile, set_ixs)
                push!(membership2tile[set_ixs], i)
            else
                membership2tile[set_ixs] = Vector{Int}([i])
            end
        end
        # build tile-to-set membership matrix
        nmasked_pertile = fill(0, length(membership2tile))
        nunmasked_pertile = fill(0, length(membership2tile))
        set2tile_ixs = [Vector{Int}() for _ in 1:length(orig_setixs)]
        @inbounds for (tile_ix, (tile_set_ixs, old_tile_ixs)) in enumerate(membership2tile)
            for set_ix in tile_set_ixs
                push!(set2tile_ixs[set_ix], tile_ix)
            end
            for old_tile_ix in old_tile_ixs
                oldtile_elms = view(mosaic.elmXtile, :, old_tile_ix)
                n_oldtile_masked = 0
                for eix in oldtile_elms
                    if elmask[eix]
                        n_oldtile_masked += 1
                    end
                end
                nmasked_pertile[tile_ix] += n_oldtile_masked
                nunmasked_pertile[tile_ix] += length(oldtile_elms) - n_oldtile_masked
            end
        end
        tileXset = SparseMaskMatrix(length(nmasked_pertile), length(orig_setixs), set2tile_ixs)
        nmasked_perset = fill(0, size(tileXset, 2))
        nunmasked_perset = Vector{Int}(size(tileXset, 2))
        @inbounds for set_ix in eachindex(nmasked_perset)
            nmasked_perset[set_ix] = sum(nmasked_pertile[view(tileXset, :, set_ix)])
            nunmasked_perset[set_ix] = mosaic.set_sizes[orig_setixs[set_ix]] - nmasked_perset[set_ix]
        end
        #println("tileXset=$tileXset")
        new{T,S}(mosaic, elmask, orig_setixs, tileXset,
                 sum(elmask), nmasked_pertile, nunmasked_pertile,
                 nmasked_perset, nunmasked_perset)
    end
end

mask(mosaic::SetMosaic, mask::Union{BitVector,Vector{Bool}}) = MaskedSetMosaic(mosaic, mask)
mask{T,S}(mosaic::SetMosaic{T,S}, sel::Set{T}) = mask(mosaic, Bool[in(e, sel) for e in mosaic.ix2elm])
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
    old2new_tile_ixs = cumsum(tile_mask)
    old_tile_ixs = find(tile_mask)
    mosaic.tileXset = SparseMaskMatrix(length(old_tile_ixs), tileXset.n,
                                       tileXset.colptr, old2new_tile_ixs[tileXset.rowval]) # remove unused tiles
    mosaic.nmasked_pertile = mosaic.nmasked_pertile[old_tile_ixs]
    mosaic.nunmasked_pertile = mosaic.nunmasked_pertile[old_tile_ixs]
    mosaic.nmasked_perset = mosaic.nmasked_perset[setmask]
    mosaic.nunmasked_perset = mosaic.nunmasked_perset[setmask]
    return mosaic
end
