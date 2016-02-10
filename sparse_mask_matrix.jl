immutable SparseMaskMatrix
    m::Int
    n::Int
    colptr::Vector{Int}
    rowval::Vector{Int}

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

function sparse_mask(m::Integer, n::Integer, rowvals_percol::Vector{Vector{Int}})
    colptr = sizehint!(Vector{Int}(), n+1)
    rowval = Vector{Int}()
    for rowvals in rowvals_percol
        push!(colptr, length(rowval)+1)
        append!(rowval, rowvals)
    end
    push!(colptr, length(rowval)+1)
    SparseMaskMatrix(m, n, colptr, rowval)
end

function sparse_mask{T}(sets, elm2ix::Dict{T, Int})
    if isempty(sets)
        return SparseMaskMatrix(length(elm2ix), 0, fill(0, 1), Vector{Int}())
    end
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

function Base.convert(::Type{SparseMaskMatrix}, mtx::Matrix{Bool})
    colptr = sizehint!(Vector{Int}(), size(mtx, 2)+1)
    rowval = Vector{Int}()
    for j in 1:size(mtx, 2)
        push!(colptr, length(rowval)+1)
        for i in 1:size(mtx, 1)
            if mtx[i, j]
                push!(rowval, i)
            end
        end
    end
    push!(colptr, length(rowval)+1)
    return SparseMaskMatrix(size(mtx)..., colptr, rowval)
end

sparse_mask(mtx::Matrix{Bool}) = convert(SparseMaskMatrix, mtx)

Base.size(mtx::SparseMaskMatrix) = (mtx.m, mtx.n)
function Base.size(mtx::SparseMaskMatrix, dim::Integer)
    (1 <= dim <= 2) || throw(ArgumentError("SparseMaskMatrix does not have dimension #$dim"))
    dim == 1 ? mtx.m : mtx.n
end
_colrange(mtx::SparseMaskMatrix, col::Integer) = mtx.colptr[col]:(mtx.colptr[col+1]-1)

Base.getindex(mtx::SparseMaskMatrix, ::Colon, col::Integer) = mtx.rowval[_colrange(mtx, col)]
function Base.getindex(mtx::SparseMaskMatrix, ::Colon, cols::Vector{Int})
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

function Base.getindex(mtx::SparseMaskMatrix, ::Colon, colmask::Union{Vector{Bool},BitVector})
    length(colmask) == size(mtx, 2) || throw(ArgumentError("Column mask length should match the number of columns"))
    mtx[:, find(colmask)]
end

Base.slice(mtx::SparseMaskMatrix, ::Colon, col::Integer) = slice(mtx.rowval, _colrange(mtx, col))
