"""
    ArrayPool{T}

Helps to maintain the pool of reusable arrays of different sizes
and reduce the burden on garbage collection.

# Type parameters
 * `T`: the type of array elements

"""
struct ArrayPool{T}
    len_pools_lock::Threads.SpinLock        # lock for thread-safe access to the pool
    len_pools::Dict{Int, Vector{Vector{T}}} # pools of vectors of different lengths

    ArrayPool{T}() where T =
        new{T}(Threads.SpinLock(), Dict{Int, Vector{Vector{T}}}())
end

"""
    acquire!(pool::ArrayPool{T}, size) where T -> Array{T}

Gets an array of a specific size from the pool.
The acquired array should be returned back to the pool using [`release!`](@ref).
The `size` could be either an integer or a tuple of integers.
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

acquire!(pool::ArrayPool, size::Tuple) =
    reshape(acquire!(pool, prod(size)), size)

"""
    release!(pool::ArrayPool{T}, arr::Array{T}) where T

Releases the array previously obtained by [`acquire!`](@ref) back into the pool.
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
