"""
Helps to maintain the pool of reusable arrays of different sizes
and reduce the burden on garbage collection.
"""
struct ArrayPool{T}
    len_pools_lock::Threads.SpinLock
    len_pools::Dict{Int, Vector{Vector{T}}} # pools of vectors of different lengths

    ArrayPool{T}() where T =
        new{T}(Threads.SpinLock(), Dict{Int, Vector{Vector{T}}}())
end

"""
Gets an array of specific size from the pool.
The returned array should be returned back to the pool using `release!()`.
"""
function acquire!(pool::ArrayPool{T}, len::Integer) where T
    lock(pool.len_pools_lock)
    len_pool = haskey(pool.len_pools, len) ?
               pool.len_pools[len] :
               get!(() -> Vector{Vector{T}}(), pool.len_pools, len)
    res = isempty(len_pool) ? Vector{T}(undef, len) : pop!(len_pool)
    unlock(pool.len_pools_lock)
    return res
end

"""
Gets an array of specific size from the pool.
The returned array should be returned back to the pool using `release!()`.
"""
acquire!(pool::ArrayPool, size::Tuple) =
    reshape(acquire!(pool, prod(size)), size)

"""
Releases an array returned by `acquire!()` back into the pool.
"""
function release!(pool::ArrayPool{T}, arr::Array{T}) where T
    len = length(arr)
    len_pool = get(pool.len_pools, len, nothing)
    if len_pool !== nothing
        #@info "release($(size(arr))) ($(length(len_pool)))"
        (length(len_pool) <= 100) || error("Overflow of $len-sized vectors pool")
        lock(pool.len_pools_lock)
        push!(len_pool, vec(arr))
        unlock(pool.len_pools_lock)
    else
        throw(DimensionMismatch("No $len-element arrays were acquired before"))
    end
    return pool
end
