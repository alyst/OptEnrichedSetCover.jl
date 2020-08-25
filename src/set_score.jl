@inline function checktail(tail::Symbol)
    if tail ∉ (:left, :right, :both)
        throw(ArgumentError("Unsupported tail specifier ($tail)"))
    end
end

"""
    logpvalue(nisect::Integer, na::Integer, nb::Integer, ntotal::Integer,
              [tail::Symbol = :right])

Log P-value for the two sets intersection.

`A` has `na` elemnts, B has `nb` elements,
they have `nisect` elements in common, and there are
`ntotal` elements in the "universe".

`tail` controls the null hypothesis:
* `:right` (default): by chance `A` and `B` would have ≥ elements in common
* `:left`: by chance `A` and `B` would have ≤ elements in common
* `:both`: by chance `A` and `B` would have either ≤ or ≥ elements in common, whichever is less probable
"""
function logpvalue(nisect::Integer, na::Integer, nb::Integer, ntotal::Integer,
                   tail::Symbol = :right)
    ((na >= 0) && (nb >= 0) && (ntotal >= 0)) ||
        throw(ArgumentError("Sets with negative number of elements"))
    ((na <= ntotal) && (nb <= ntotal)) ||
        throw(ArgumentError("Sets bigger than total number of elements"))
    # corner cases
    isect_max = min(na, nb)
    isect_min = max(0, na + nb - ntotal)
    if tail == :right
        (nisect > isect_max) && return -Inf
        (nisect <= isect_min) && return 0.0
    elseif tail == :left
        (nisect >= isect_max) && return 0.0
        (nisect < isect_min) && return -Inf
    elseif tail != :both
        # FIXME corner cases when tail == :both
        throw(ArgumentError("Unsupported tail specifier ($tail)"))
    end
    # normal cases
    distr = Distributions.Hypergeometric(na, ntotal - na, nb)
    if tail == :right
        return nisect < isect_max ? logccdf(distr, nisect-1) : logpdf(distr, nisect)
    elseif tail == :left
        return nisect > isect_min ? logcdf(distr, nisect) : logpdf(distr, nisect)
    elseif tail == :both
        return log(2.0) + min(logcdf(distr, nisect), logccdf(distr, nisect-1), log(0.5))
    else
        return NaN # noop
    end
end

# size of the sets intersection
# the sets are represented by the sorted sets of their tiles
function _isect_size(set1_tiles::AbstractVector{Int},
                     set2_tiles::AbstractVector{Int},
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

# returns set×set score matrix
# a×b score is the Fisher's test P-value of their intersection
# a×a self-intersections are processed in the same way
# (might be used later in cover problem construction to penalize the use of
# multiple variables referencing the same set)
function _setXset_scores(tileXset::SparseMaskMatrix, total_size::Int,
                         set_sizes::AbstractVector{Int},
                         tile_sizes::AbstractVector{Int},
                         used_sets = Base.OneTo(size(tileXset, 2)))
    nsets = length(used_sets)
    res = zeros(Float64, nsets, nsets)
    # fill upper triangle
    @inbounds for i in 1:nsets
        set1_ix = used_sets[i]
        set1_size = set_sizes[set1_ix]
        set1_size == 0 && continue # skip empty sets
        set1_tiles = view(tileXset, :, set1_ix)
        for j in i:nsets
            set2_ix = used_sets[j]
            set2_size = set_sizes[set2_ix]
            set2_size == 0 && continue # skip empty sets
            set2_tiles = view(tileXset, :, set2_ix)
            # one-sided Fisher's P-value, right tail
            isect_size = _isect_size(set1_tiles, set2_tiles, tile_sizes)
            res[i, j] = min(logpvalue(isect_size, set1_size, set1_size + set2_size - isect_size, total_size),
                            logpvalue(isect_size, set2_size, set1_size + set2_size - isect_size, total_size))
        end
    end
    # symmetrize: copy upper triangle to lower
    @inbounds for i in 1:nsets
        for j in (i+1):nsets
            res[j, i] = res[i, j]
        end
    end
    return res
end

"""
    set_relevance(nset_observed::Integer, nset::Integer,
                  nobserved::Integer, ntotal::Integer) -> Float64

Calculates the relevance weight of the set that contains `nset` elements,
`nset_observed` of which were present (not necessarily enriched) in the data
that identified `nobserved` elements out of all known (`ntotal`).
It is used by [`SetMosaic`](@ref) to penalize the sets,
which could not be observed in the data (e.g. biological processes or pathways
that involve proteins not expressed by the cells used in the experiments).

While for [`MaskedSetMosaic`](@ref) it's recommended to use the IDs of data entities
(e.g. protein group IDs for proteomic data) to correctly count the set sizes
and estimate enrichment; `set_relevance()` should use the counts derived from
the original IDs of the annotation database (e.g. UniProt accession codes).
Otherwise it's not possible to correctly estimate the number of elements that
belong to the given annotated set, but were not observed in the data.

The returned value is the probability that no more than `nset_observed`
elements were observed at random.
"""
set_relevance(nset_observed::Integer, nset::Integer,
              nobserved::Integer, ntotal::Integer) =
    cdf(Hypergeometric(nset, ntotal - nset, nobserved), nset_observed)
