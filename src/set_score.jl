@inline function checktail(tail::Symbol)
    if tail ∉ (:left, :right, :both)
        throw(ArgumentError("Unsupported tail specifier ($tail)"))
    end
end

"""
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
    ((na >= 0) && (nb >= 0) && (ntotal >= 0)) || throw(ArgumentError("Sets with negative number of elements"))
    ((na <= ntotal) && (nb <= ntotal)) || throw(ArgumentError("Sets bigger than total number of elements"))
    # corner cases
    if nisect > min(na, nb)
        checktail(tail)
        return tail == :left ? 0.0 : -Inf
    elseif nisect < min(0, na + nb - ntotal)
        checktail(tail)
        return tail == :right ? 0.0 : -Inf
    elseif ntotal==max(na,nb) && nisect==min(na,nb) # FIXME remove this corner case when Rmath-julia (and R upstream would be fixed)
        checktail(tail)
        return 0.0
    end
    # normal cases
    distr = Distributions.Hypergeometric(na, ntotal - na, nb)
    if tail == :right
        return logccdf(distr, nisect-1)
    elseif tail == :left
        return logcdf(distr, nisect)
    elseif tail == :both
        return log(2.0) + min(logcdf(distr, nisect), logccdf(distr, nisect-1), log(0.5))
    else
        checktail(tail)
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

function _setXset_scores(tileXset::SparseMaskMatrix, total_size::Int,
                         set_sizes::AbstractVector{Int},
                         tile_sizes::AbstractVector{Int},
                         used_sets = Base.OneTo(size(tileXset, 2)))
    nsets = length(used_sets)
    res = zeros(Float64, nsets, nsets)
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
            res[i, j] =
                logpvalue(_isect_size(set1_tiles, set2_tiles, tile_sizes),
                          set1_size, set2_size, total_size)
        end
    end
    # symmetrize
    @inbounds for i in 1:nsets
        for j in (i+1):nsets
            res[j, i] = res[i, j]
        end
    end
    return res
end
