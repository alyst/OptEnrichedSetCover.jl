"""
Sparse representation of `Matrix{Bool}`, when all "falses" are "structural".
It uses compressed sparse column (CSC) representation, except no
non-zero values need to be stored.
"""
struct SparseMaskMatrix
    m::Int              # rows number
    n::Int              # cols number
    colptr::Vector{Int} # indices of column starts in rowval vector
    rowval::Vector{Int} # indices of rows of true matrix elements

    function SparseMaskMatrix(m::Integer, n::Integer, colptr::Vector{Int}, rowval::Vector{Int})
        length(colptr) == n+1 || throw(ArgumentError("col pointers ($(length(colptr))) does not match columns count ($n)"))
        colptr[end] == length(rowval)+1 || throw(ArgumentError("colptr[end] does not match rowval length"))
        new(Int(m), Int(n), colptr, rowval)
    end

    function SparseMaskMatrix(m::Integer=0, n::Integer=0) # empty mask
        new(m, n, fill(0, n+1), Vector{Int}())
    end
end

Base.copy(mtx::SparseMaskMatrix) = SparseMaskMatrix(mtx.m, mtx.n, copy(mtx.colptr), copy(mtx.rowval))

"""
Construct `SparseMaskMatrix` given vector of true row indices per each column.
"""
function SparseMaskMatrix(m::Integer, rowvals_percol::Vector{Vector{Int}})
    colptr = sizehint!(Vector{Int}(), length(rowvals_percol)+1)
    nrowvals_cumlenp1 = 1
    for rowvals in rowvals_percol
        push!(colptr, nrowvals_cumlenp1)
        nrowvals_cumlenp1 += length(rowvals)
    end
    push!(colptr, nrowvals_cumlenp1)
    allrowvals = sizehint!(Vector{Int}(), nrowvals_cumlenp1-1)
    for rowvals in rowvals_percol
        append!(allrowvals, rowvals)
    end
    SparseMaskMatrix(m, length(rowvals_percol), colptr, allrowvals)
end

"""
Construct `SparseMaskMatrix` from the family of sets.

 * `sets` family of sets, one set per result column
 * `elm2ix` mapping from set element to its index (mask row index)
"""
function SparseMaskMatrix(sets, elm2ix::Dict{T, Int}) where T
    isempty(sets) && return SparseMaskMatrix(length(elm2ix), 0, fill(0, 1), Vector{Int}())

    elm_ixs = Vector{Int}()
    set_ranges = Vector{Int}()
    sizehint!(set_ranges, length(sets))
    for s in sets
        push!(set_ranges, length(elm_ixs)+1)
        append!(elm_ixs, sort!([elm2ix[e] for e in s]))
    end
    push!(set_ranges, length(elm_ixs)+1)
    return SparseMaskMatrix(length(elm2ix), length(sets), set_ranges, elm_ixs)
end

function SparseMaskMatrix(mtx::AbstractMatrix{Bool})
    colptr = sizehint!(Vector{Int}(), size(mtx, 2)+1)
    rowval = Vector{Int}()
    @inbounds for j in 1:size(mtx, 2)
        push!(colptr, length(rowval)+1)
        for i in 1:size(mtx, 1)
            mtx[i, j] && push!(rowval, i)
        end
    end
    push!(colptr, length(rowval)+1)
    return SparseMaskMatrix(size(mtx)..., colptr, rowval)
end

Base.size(mtx::SparseMaskMatrix) = (mtx.m, mtx.n)
function Base.size(mtx::SparseMaskMatrix, dim::Integer)
    (1 <= dim <= 2) || throw(ArgumentError("SparseMaskMatrix does not have dimension #$dim"))
    dim == 1 ? mtx.m : mtx.n
end
@inline _colrange(mtx::SparseMaskMatrix, col::Integer) = mtx.colptr[col]:(mtx.colptr[col+1]-1)

@inline Base.getindex(mtx::SparseMaskMatrix, ::Colon, col::Integer) = mtx.rowval[_colrange(mtx, col)]

Base.@propagate_inbounds function Base.getindex(mtx::SparseMaskMatrix, ::Colon, cols::AbstractVector{Int})
    rowvals = Vector{Int}()
    colptr = Vector{Int}()
    for col in cols
        # FIXME avoid temp array alloc
        push!(colptr, length(rowvals)+1)
        append!(rowvals, mtx[:, col])
    end
    push!(colptr, length(rowvals)+1)
    SparseMaskMatrix(mtx.m, length(cols), colptr, rowvals)
end

Base.@propagate_inbounds function Base.getindex(mtx::SparseMaskMatrix, rows::AbstractVector{Int}, ::Colon)
    old2new = Dict(oldix => newix for (newix, oldix) in enumerate(rows))
    rowval = Vector{Int}()
    thiscol = Vector{Int}()
    colptr = [1]
    @inbounds for j in 1:mtx.n
        # collect new row indices of j-th column
        empty!(thiscol)
        for i in _colrange(mtx, j)
            newrow = get(old2new, mtx.rowval[i], 0)
            (newrow != 0) && push!(thiscol, newrow)
        end
        append!(rowval, sort!(thiscol))
        push!(colptr, length(rowval)+1)
    end
    SparseMaskMatrix(length(rows), mtx.n, colptr, rowval)
end

Base.@propagate_inbounds function Base.getindex(mtx::SparseMaskMatrix, ::Colon, colmask::AbstractVector{Bool})
    length(colmask) == size(mtx, 2) || throw(ArgumentError("Column mask length should match the number of columns"))
    mtx[:, find(colmask)]
end

Base.view(mtx::SparseMaskMatrix, ::Colon, col::Integer) = view(mtx.rowval, _colrange(mtx, col))

function Base.convert(::Type{Matrix{Bool}}, mtx::SparseMaskMatrix)
    res = fill(false, size(mtx))
    @inbounds for i in 1:size(mtx, 2)
        for j in _colrange(mtx, i)
            res[mtx.rowval[j], i] = true
        end
    end
    return res
end

Base.convert(::Type{Matrix}, mtx::SparseMaskMatrix) = convert(Matrix{Bool}, mtx)
