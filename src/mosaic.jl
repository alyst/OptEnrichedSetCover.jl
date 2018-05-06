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
            ix = get(elm2ix, e, 0)
            if ix > 0
                @inbounds set_elms[ix] = true
            end
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

# unlike Base.union() does not use recursion
function _union{T}(::Type{T}, sets)
    res = Set{T}()
    for set in sets
        union!(res, set)
    end
    return res
end

function _set_sizes(tileXset::SparseMaskMatrix, tile_sizes::Vector{Int})
    set_sizes = zeros(Int, size(tileXset, 2))
    for six in eachindex(set_sizes)
        tile_ixs = view(tileXset, :, six)
        isempty(tile_ixs) && continue
        @inbounds set_sizes[six] = sum(tix -> tile_sizes[tix], view(tileXset, :, six))
    end
    return set_sizes
end

"""
A collection of (potentially overlapping) sets as
a "mosaic" of non-overlapping "tiles".

* `T` type of elements
* `S` type of set keys
"""
struct SetMosaic{T,S}
    ix2elm::Vector{T}           # element index to element
    elm2ix::Dict{T, Int}        # element to its index

    ix2set::Vector{S}           # set index to the ID of the most relevant set (if there are duplicates)
    set2ix::Dict{S, Int}        # set ID to index

    set_sizes::Vector{Int}
    set_relevances::Vector{Float64}  # 0..1 set relevance score

    setXelm::SparseMaskMatrix   # set×element membership
    elmXset::SparseMaskMatrix   # element×set membership (transpose of setXelm)
    elmXtile::SparseMaskMatrix  # element-to-tile mask
    tileXset::SparseMaskMatrix  # rows = tiles, cols = sets from collection, true if tile is a subset of a set

    setXset_scores::Matrix{Float64}

    function SetMosaic(setids::AbstractVector{S},
                       sets::AbstractVector{<:Union{AbstractVector{T},Set{T}}},
                       all_elms::Set{T},
                       set_relevances::AbstractVector{<:Real}) where {T, S}
        length(setids) == length(sets) ||
            throw(ArgumentError("Number of set IDs does not match the number of sets"))
        length(set_relevances) == length(sets) ||
            throw(ArgumentError("Number of set relevance scores does not match the number of sets"))
        ix2elm, elm2ix = _encode_elements(all_elms)
        setXelm, elmXset, elmXtile, tileXset = _prepare_tiles(sets, elm2ix)

        # elementwise sets comparison
        function setels_compare(iset, jset)
            lendiff = length(iset) - length(jset) # prioritize length since we need sorting just to find duplicated
            (lendiff != 0) && return lendiff > 0 ? 1 : -1
            (length(iset) == 0) && return 0 # both empty
            for i in eachindex(iset)
                @inbounds eldiff = iset[i] - jset[i]
                (eldiff != 0) && return eldiff > 0 ? 1 : -1
            end
            return 0
        end
        # sort sets by their tiles and relevance, so that duplicate sets are adjacent
        function sets_isless(i, j)
            els_diff = setels_compare(view(tileXset, :, i), view(tileXset, :, j))
            (els_diff != 0) && return els_diff < 0
            rel_diff = set_relevances[i] - set_relevances[j]
            (rel_diff != 0.0) && return rel_diff > 0.0 # descending relevance
            return setids[i] < setids[j]
        end
        @inbounds set_sort = sortperm(1:size(tileXset, 2), lt=sets_isless)
        ix2new = fill(0, length(setids))
        new2ix = Vector{Int}()
        if !isempty(set_sort) # remove duplicate sets keeping the ones with higher relevance and id
            lastset_tiles = view(tileXset, :, 1)
            for setix in set_sort
                curset_tiles = view(tileXset, :, setix)
                if isempty(new2ix) || curset_tiles != lastset_tiles
                    lastset_tiles = curset_tiles
                    push!(new2ix, setix)
                end
                ix2new[setix] = length(new2ix)
            end
            neword = sortperm(new2ix)
            new2ix = new2ix[neword]
            ix2new = invperm(neword)[ix2new]
        end

        tile_sizes = Int[length(view(elmXtile, :, tile_ix)) for tile_ix in 1:size(elmXtile, 2)]
        tileXset = tileXset[:, new2ix]
        set_sizes = _set_sizes(tileXset, tile_sizes)
        new{T, S}(ix2elm, elm2ix, setids[new2ix],
                  Dict(id => ix2new[ix] for (ix, id) in enumerate(setids)),
                  set_sizes, set_relevances[new2ix],
                  setXelm[new2ix, :], elmXset[:, new2ix], elmXtile, tileXset,
                  _setXset_scores(tileXset, length(all_elms), set_sizes, tile_sizes))
    end
end

"""
Construct `SetMosaic` for a given nameless sets collection.
"""
SetMosaic(sets::AbstractVector{Set{T}}, all_elms::Set{T} = _union(T, sets),
          set_relevances::AbstractVector{Float64} = fill(1.0, length(sets))) where {T} =
    SetMosaic(collect(1:length(sets)), sets, all_elms, set_relevances)

"""
Constructs `SetMosaic` for a given named sets collection.
"""
SetMosaic(sets::Dict{S, Set{T}}, all_elms::Set{T} = _union(T, values(sets)),
          set_relevances::Dict{S, Float64} = Dict(k => 1.0 for k in keys(sets))) where {S, T} =
    SetMosaic(collect(keys(sets)), collect(values(sets)), all_elms,
              [set_relevances[k] for k in keys(sets)])

nelements(mosaic::SetMosaic) = length(mosaic.ix2elm)
ntiles(mosaic::SetMosaic) = size(mosaic.tileXset, 1)
nsets(mosaic::SetMosaic) = size(mosaic.tileXset, 2)

tile(mosaic::SetMosaic, tile_ix::Integer) = view(mosaic.elmXtile, :, tile_ix)
set(mosaic::SetMosaic, set_ix::Integer) = view(mosaic.elmXset, :, set_ix)

setsize(mosaic::SetMosaic, set_ix::Integer) = mosaic.set_sizes[set_ix]
