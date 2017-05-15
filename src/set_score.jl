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
function logpvalue(na::Integer, nb::Integer, ntotal::Integer, nisect::Integer,
                   tail::Symbol = :right)
    ((na >= 0) && (nb >= 0) && (ntotal >= 0)) || throw(ArgumentError("Sets with negative number of elements"))
    ((na <= ntotal) && (nb <= ntotal)) || throw(ArgumentError("Sets bigger than total number of elements"))
    # corner cases
    if nisect > min(na, nb)
        return tail == :left ? 0.0 : -Inf
    elseif nisect < min(0, na + nb - ntotal)
        return tail == :right ? 0.0 : -Inf
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
        throw(ArgumentError("Unsupported tail specifier ($tail)"))
    end
end
