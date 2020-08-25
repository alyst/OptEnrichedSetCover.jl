# assign index (1..n) to each elements of the set
# returns the tuple:
#   - the vector of elements ordered by their index
#   - the dictionary from element to its index
function _index_elements(all_elms::Set{T}) where T
    ix2elm = collect(all_elms)
    sort!(ix2elm)
    # map from element to its index
    elm2ix = Dict{T,Int}()
    sizehint!(elm2ix, length(ix2elm))
    for i in 1:length(ix2elm)
        elm2ix[ix2elm[i]] = i
    end
    return ix2elm, elm2ix
end

# prepares disjoint tiles for the sets collection
# returns the tuple of :
#   - transposed version of the next sparse mask (set×element)
#   - sparse mask representation of sets collection (element×set)
#   - sparse mask representation of tiles collecttion (element×tile)
#   - sparse mask representation of sets composition using tiles (tile×set)
function _prepare_tiles(sets, elm2ix::Dict{<:Any, Int})
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
        @inbounds set_ixs = findall(view(setXelm, :, i))
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

    return SparseMaskMatrix(setXelm), SparseMaskMatrix(elmXset),
        SparseMaskMatrix(length(elm2ix), length(tile_elm_ranges)-1,
                         tile_elm_ranges, elm_ixs), # elmXtile
        SparseMaskMatrix(length(tile_elm_ranges)-1, length(sets),
                         set_tile_ranges, tile_ixs) # tileXset
end

# calculate the sizes of sets given the sizes of tiles
function _set_sizes(tileXset::SparseMaskMatrix,
                    tile_sizes::AbstractVector{Int})
    @assert length(tile_sizes) == size(tileXset, 1)
    set_sizes = zeros(Int, size(tileXset, 2))
    for six in eachindex(set_sizes)
        tile_ixs = view(tileXset, :, six)
        isempty(tile_ixs) && continue
        @inbounds set_sizes[six] = sum(tix -> tile_sizes[tix], tile_ixs)
    end
    return set_sizes
end

"""
    SetMosaic{T,S}

Represents a collection of (potentially overlapping) sets as
a "mosaic" of non-overlapping "tiles".

# Type parameters
* `T`: type of set elements
* `S`: type of set keys
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
                       set_relevances::Union{AbstractVector{<:Real}, Nothing} = nothing;
                       setXset_nextra_elms::Integer=0) where {T, S}
        length(setids) == length(sets) ||
            throw(ArgumentError("Number of set IDs does not match the number of sets"))
        (set_relevances === nothing) || (length(set_relevances) == length(sets)) ||
            throw(ArgumentError("Number of set relevance scores does not match the number of sets"))
        ix2elm, elm2ix = _index_elements(all_elms)
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
            rel_diff = set_relevances === nothing ? 0.0 : set_relevances[i] - set_relevances[j]
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
                  set_sizes, set_relevances === nothing ? fill(1.0, length(new2ix)) : set_relevances[new2ix],
                  setXelm[new2ix, :], elmXset[:, new2ix], elmXtile, tileXset,
                  _setXset_scores(tileXset, length(all_elms) + setXset_nextra_elms, set_sizes, tile_sizes))
    end
end

"""
    SetMosaic(sets, [all_elms::Set{T}],
              [set_relevances::AbstractVector{Float64}];
              [setXset_nextra_elms=0]) where {T, S} -> SetMosaic{T, S}

Construct `SetMosaic` for a given sets collection.

# Arguments
  * `sets`: a mapping from the set ids (of type `S`) to its elements (`Set{T}`).
     Could be either a dictionary or a vector (then *S* is *Int*).
  * `setXset_nextra_elms`: added to the number of total elements to adjust
    redundancy scores (so that overlaps of very big terms are still penalized)
"""
SetMosaic(sets::AbstractVector{Set{T}}, all_elms::AbstractSet{T} = foldl(union!, sets, init=Set{T}()),
          set_relevances::Union{AbstractVector{Float64}, Nothing} = nothing;
          kwargs...) where {T} =
    SetMosaic(collect(1:length(sets)), sets, all_elms, set_relevances; kwargs...)

SetMosaic(sets::Dict{S, Set{T}}, all_elms::Set{T} = foldl(union!, values(sets), init=Set{T}()),
          set_relevances::Union{Dict{S, Float64}, Nothing} = nothing;
          kwargs...) where {S, T} =
    SetMosaic(collect(keys(sets)), collect(values(sets)), all_elms,
              set_relevances === nothing ? nothing :
              get.(Ref(set_relevances), keys(sets), 1.0); kwargs...)

"""
    nelements(mosaic::SetMosaic) -> Int

The number of distinct elements (i.e. genes) in the sets collection.
"""
nelements(mosaic::SetMosaic) = length(mosaic.ix2elm)

"""
    nsets(mosaic::SetMosaic) -> Int

The number of tiles (pairwise disjoint sets) in the mosaic representation of the set collection.
"""
ntiles(mosaic::SetMosaic) = size(mosaic.tileXset, 1)

"""
    nsets(mosaic::SetMosaic) -> Int

The number of sets in the collection.
"""
nsets(mosaic::SetMosaic) = size(mosaic.tileXset, 2)

"""
    tile(mosaic::SetMosaic, i::Integer) -> AbstractVector{Int}

Get the `i`-th mosaic tile as the vector of indices of its elements.
"""
tile(mosaic::SetMosaic, tile_ix::Integer) = view(mosaic.elmXtile, :, tile_ix)

"""
    set(mosaic::SetMosaic, i::Integer) -> AbstractVector{Int}

Get the `i`-th set as the vector of indices of its elements.
"""
set(mosaic::SetMosaic, set_ix::Integer) = view(mosaic.elmXset, :, set_ix)

"""
    set(mosaic::SetMosaic, i::Integer) -> Int

Get the size of the `i`-th set.
"""
setsize(mosaic::SetMosaic, set_ix::Integer) = mosaic.set_sizes[set_ix]
