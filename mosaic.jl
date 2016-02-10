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
        return Matrix{Bool}(0,0), Vector{Vector{Int}}()
    end

    # build set-X-element membership matrix
    setXelm = fill(false, length(sets), length(elm2ix))
    for (i, s) in enumerate(sets)
        for e in s
            setXelm[i, elm2ix[e]] = true
        end
    end
    # squash element-X-set matrix (duplicate rows) into membership-mask-to-tile dict
    membership2tile = Dict{Vector{Bool}, Vector{Int}}()
    sizehint!(membership2tile, length(elm2ix))
    for i in 1:size(setXelm, 2)
        mask = setXelm[:, i]
        if haskey(membership2tile, mask)
            push!(membership2tile[mask], i)
        else
            membership2tile[mask] = fill(i, 1)
        end
    end
    # build tile-X-set membership matrix
    tileXset = fill(false, length(membership2tile), length(sets))
    tiles = Vector{Vector{Int}}()
    sizehint!(tiles, length(tileXset))
    for (mask, tile) in membership2tile
        push!(tiles, tile)
        tileXset[length(tiles), :] = mask
    end
    tileXset, tiles
end

function _setXset_scores(tileXset::Matrix{Bool}, total_size::Int, set_sizes::Vector{Int}, tile_sizes::Vector{Int})
    res = Matrix{Float64}(size(tileXset, 2), size(tileXset, 2))
    @inbounds for set1_ix in 1:size(tileXset, 2)
        res[set1_ix, set1_ix] = 0.0
        for set2_ix in (set1_ix+1):size(tileXset, 2)
            isect_size = 0
            for tile_ix in eachindex(tile_sizes)
                if tileXset[tile_ix, set1_ix] && tileXset[tile_ix, set2_ix]
                    isect_size += tile_sizes[tile_ix]
                end
            end
            # one-sided Fisher's P-value
            res[set1_ix, set2_ix] = res[set2_ix, set1_ix] = logpvalue(set_sizes[set1_ix], set_sizes[set2_ix], total_size, isect_size, tail=:left);
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
  Encodes a collection of (overlapping) sets as
  a mosaic of non-overlapping "tiles".

  * `T` type of elements
  * `S` type of set keys
"""
immutable SetMosaic{T,S}
    ix2elm::Vector{T}       # element index to element
    elm2ix::Dict{T, Int}    # element to its index

    ix2set::Vector{S}       # set index to ID
    set2ix::Dict{S, Int}    # set ID to index

    set_sizes::Vector{Int}

    tileXset::Matrix{Bool}  # rows = tiles, cols = sets from collection, true if tile is a subset of a set
    tiles::Vector{Vector{Int}} # tile elements

    setXset_scores::Matrix{Float64}

    """
      Constructs `SetMosaic` for a given nameless sets collection.
    """
    function Base.call{T}(::Type{SetMosaic}, sets::Vector{Set{T}}, all_elms::Set{T} = _union(T, sets))
        ix2elm, elm2ix = _encode_elements(all_elms)
        tileXset, tiles = _prepare_tiles(sets, elm2ix)
        tile_sizes = Int[length(tile) for tile in tiles]
        set_sizes = Int[length(set) for set in sets]
        new{T, Int}(ix2elm, elm2ix,
                    collect(eachindex(sets)), Dict(Pair{Int,Int}[Pair(i,i) for i in eachindex(sets)]),
                    set_sizes,
                    tileXset, tiles, _setXset_scores(tileXset, length(all_elms), set_sizes, tile_sizes))
    end

    """
      Constructs `SetMosaic` for a given named sets collection.
    """
    function Base.call{T, S}(::Type{SetMosaic}, sets::Dict{S, Set{T}}, all_elms::Set{T} = _union(T, values(sets)))
        ix2elm, elm2ix = _encode_elements(all_elms)
        tileXset, tiles = _prepare_tiles(values(sets), elm2ix)
        tile_sizes = Int[length(tile) for tile in tiles]
        set_sizes = Int[length(set) for set in values(sets)]
        new{T, S}(ix2elm, elm2ix,
                  collect(keys(sets)), Dict{S, Int}(Pair{S,Int}[Pair(s, i) for (i, s) in enumerate(keys(sets))]),
                  set_sizes,
                  tileXset, tiles, _setXset_scores(tileXset, length(all_elms), set_sizes, tile_sizes))
    end
end

nelements(mosaic::SetMosaic) = length(mosaic.ix2elm)
ntiles(mosaic::SetMosaic) = length(mosaic.tiles)
nsets(mosaic::SetMosaic) = size(mosaic.tileXset, 2)
tiles(mosaic::SetMosaic) = mosaic.tiles

"""
    `SetMosaic` with an elements mask on top.
    Sets that are not overlapping with the mask are excluded from `MaskedSetMosaic`.
"""
type MaskedSetMosaic{T,S}
    original::SetMosaic{T,S}    # original mosaic
    mask::Vector{Bool}          # elements mask
    setixs::Vector{Int}         # original indices of the included sets

    """
      Subset of the original mosaic that excludes sets not overlapping with
      the mask and their tiles.
      The remaining tiles that have the same membership within the remaining

      are squashed.
    """
    tileXset::Matrix{Bool}
    tiles::Vector{Vector{Int}} # tile elements
    total_masked::Int      # total number of masked elements
    nmasked_pertile::Vector{Int}   # number of masked elements in a tile
    nunmasked_pertile::Vector{Int} # number of unmasked elements in a tile
    nmasked_perset::Vector{Int}   # number of masked elements in a set
    nunmasked_perset::Vector{Int} # number of unmasked elements in a set

    function MaskedSetMosaic(original::SetMosaic{T, S}, mask::Vector{Bool},
                             setixs::Vector{Int},
                             tileXset::Matrix{Bool}, tiles::Vector{Vector{Int}},
                             total_masked::Number,
                             nmasked_pertile::Vector{Int}, nunmasked_pertile::Vector{Int},
                             nmasked_perset::Vector{Int}, nunmasked_perset::Vector{Int}
    )
        new(original, mask, setixs, tileXset, tiles, total_masked,
            nmasked_pertile, nunmasked_pertile,
            nmasked_perset, nunmasked_perset)
    end

    function Base.call{T,S}(::Type{MaskedSetMosaic}, mosaic::SetMosaic{T, S}, mask::Vector{Bool})
        length(mask) == nelements(mosaic) ||
            throw(ArgumentError("Mask length ($(length(mask))) should match the number of elements ($(nelements(mosaic)))"))
        #println("mask=", mask)
        # get the tiles that overlap with the mask
        # with tiles it's faster than with the sets because no element is repeated
        tileixs = Vector{Int}()
        for i in 1:length(mosaic.tiles)
            #println("Tile #$i")
            for e in mosaic.tiles[i]
                if mask[e]
                    #println("Masked element $e is in tile $i: $(mosaic.tiles[i])")
                    push!(tileixs, i)
                    break
                end
            end
        end
        #println("masked tiles=$tileixs")
        # mask of sets that overlap with the masked elements
        setmask = squeeze(any(mosaic.tileXset[tileixs, :], 1), 1)
        #println("setmask=$setmask")
        subsetXtile = mosaic.tileXset[:, setmask]'
        #println("subsetXtile=$subsetXtile")
        # detect and squash tiles that have the same membership pattern wrt the masked sets
        membership2tile = Dict{Vector{Bool}, Vector{Int}}()
        sizehint!(membership2tile, size(subsetXtile, 2))
        for i in 1:size(subsetXtile, 2)
            tile_mask = subsetXtile[:, i]
            if haskey(membership2tile, tile_mask)
                push!(membership2tile[tile_mask], i)
            else
                membership2tile[tile_mask] = Vector{Int}([i])
            end
        end
        # build tile-to-set membership matrix
        tileXset = fill(false, length(membership2tile), size(subsetXtile, 1))
        tiles = Vector{Vector{Int}}()
        nmasked_pertile = Vector{Int}(length(membership2tile))
        nunmasked_pertile = Vector{Int}(length(membership2tile))
        sizehint!(tiles, length(tileXset))
        for (tile_mask, tile_ixs) in membership2tile
            tile_ix = length(tiles)+1
            push!(tiles, sort!(vcat(mosaic.tiles[tile_ixs]...)))
            tileXset[tile_ix, :] = tile_mask
            nmasked_pertile[tile_ix] = sum(mask[tiles[tile_ix]])
            nunmasked_pertile[tile_ix] = length(tiles[tile_ix]) - nmasked_pertile[tile_ix]
        end
        nmasked_perset = fill(0, size(tileXset, 2))
        nunmasked_perset = Vector{Int}(size(tileXset, 2))
        setixs = find(setmask)
        for set_ix in eachindex(nmasked_perset)
            nmasked_perset[set_ix] = sum(nmasked_pertile[tileXset[:,set_ix]])
            nunmasked_perset[set_ix] = mosaic.set_sizes[setixs[set_ix]] - nmasked_perset[set_ix]
        end
        #println("tileXset=$tileXset")
        new{T,S}(mosaic, mask, setixs, tileXset, tiles,
                 sum(mask), nmasked_pertile, nunmasked_pertile,
                 nmasked_perset, nunmasked_perset)
    end
end

mask(mosaic::SetMosaic, mask::Vector{Bool}) = MaskedSetMosaic(mosaic, mask)
mask{T,S}(mosaic::SetMosaic{T,S}, sel::Set{T}) = mask(mosaic, Bool[in(e, sel) for e in mosaic.ix2elm])
unmask(mosaic::MaskedSetMosaic) = mosaic.original

nelements(mosaic::MaskedSetMosaic) = nelements(mosaic.original)
tiles(mosaic::MaskedSetMosaic) = mosaic.tiles
ntiles(mosaic::MaskedSetMosaic) = length(mosaic.tiles)
nsets(mosaic::MaskedSetMosaic) = size(mosaic.tileXset, 2)
nmasked(mosaic::MaskedSetMosaic) = mosaic.total_masked
nunmasked(mosaic::MaskedSetMosaic) = nelements(mosaic) - mosaic.total_masked
nmasked_pertile(mosaic::MaskedSetMosaic) = mosaic.nmasked_pertile
nunmasked_pertile(mosaic::MaskedSetMosaic) = mosaic.nunmasked_pertile
nmasked_perset(mosaic::MaskedSetMosaic) = mosaic.nmasked_perset
nunmasked_perset(mosaic::MaskedSetMosaic) = mosaic.nunmasked_perset

function Base.copy{T,S}(mosaic::MaskedSetMosaic{T,S})
    # copy everything, except the original mosaic (leave the reference to the same object)
    MaskedSetMosaic{T,S}(mosaic.original, copy(mosaic.mask), copy(mosaic.setixs),
                    copy(mosaic.tileXset), copy(mosaic.tiles), mosaic.total_masked,
                    copy(mosaic.nmasked_pertile), copy(mosaic.nunmasked_pertile),
                    copy(mosaic.nmasked_perset), copy(mosaic.nunmasked_perset))
end

"""
  Excludes sets from the mosaic and updates the set of active tiles.

  * `setmask` mask of sets to exclude
"""
function Base.filter!(mosaic::MaskedSetMosaic, setmask::Union{Vector{Bool},BitVector})
    nsets(mosaic) == length(setmask) ||
        throw(ArgumentError("Mask length ($(lenght(setmask))) does not match the number of sets in mosaic ($(nsets(mosaic)))"))
    mosaic.setixs = mosaic.setixs[setmask]
    mosaic.tileXset = mosaic.tileXset[:, setmask]
    tile_ixs = find(squeeze(any(mosaic.tileXset, 2), 2))
    mosaic.tileXset = mosaic.tileXset[tile_ixs, :] # remove unused tiles
    mosaic.tiles = mosaic.tiles[tile_ixs]
    mosaic.nmasked_pertile = mosaic.nmasked_pertile[tile_ixs]
    mosaic.nunmasked_pertile = mosaic.nunmasked_pertile[tile_ixs]
    mosaic.nmasked_perset = mosaic.nmasked_perset[setmask]
    mosaic.nunmasked_perset = mosaic.nunmasked_perset[setmask]
    return mosaic
end
