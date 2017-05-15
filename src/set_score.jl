"""
Log P-value for the two sets intersection.

`A` has `a_size` elemnts, B has `b_size` elements,
they have `isect_size` elements in common, and there are
`all_size` elements in the "universe" in total.

`tail` controls the null hypothesis:
* `:left`: by chance `A` and `B` would have more elements in common
* `:right`: by chance `A` and `B` would have less elements in common
* `:both`: by chance `A` and `B` would have either more or less elements in common
"""
function logpvalue(a_size::Integer, b_size::Integer,
                   all_size::Integer, isect_size::Integer,
                   tail::Symbol = :left)
    ((a_size >= 0) && (b_size >= 0) && (all_size >= 0)) || throw(ArgumentError("Sets with negative number of elements"))
    ((a_size <= all_size) && (b_size <= all_size)) || throw(ArgumentError("Sets bigger than total number of elements"))
    # corner cases
    if isect_size >= min(a_size, b_size)
        return tail == :right ? -Inf : 0.0
    elseif isect_size < min(0, a_size + b_size - all_size)
        return tail == :left ? -Inf : 0.0
    end
    # normal cases
    distr = Distributions.Hypergeometric(a_size, all_size - a_size, b_size)
    if tail == :left
        return logccdf(distr, isect_size-1)
    elseif tail == :right
        return logcdf(distr, isect_size)
    elseif tail == :both
        return 2.0 * min(logcdf(distr, isect_size),
                         logccdf(distr, isect_size-1), 0.5)
    else
        throw(ArgumentError("Unsupported tail specifier ($tail)"))
    end
end
