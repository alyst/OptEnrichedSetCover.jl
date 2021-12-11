# \sum_{i, j} A_{i, j} min(w_i, w_j)
function minplus_quad(A::AbstractMatrix{T},
                      w::AbstractVector{T}) where T
    length(w) == size(A, 1) == size(A, 2) ||
        throw(DimensionMismatch("A $(size(A)) and w ($(length(w))) size mismatch"))
    isempty(w) && return zero(T)
    res = zero(T)
    @inbounds for i in eachindex(w)
        wi = w[i]
        ioffset = (i-1)*size(A, 1)
        resi = zero(T)
        if zero(T) < wi < one(T)
            for j in eachindex(w)
                resi += A[ioffset+j]*min(wi, w[j])
            end
        elseif wi == one(T) # wi >= maxw
            for j in eachindex(w)
                resi += A[ioffset+j] * w[j] # FIXME dotprod() faster?
            end
        end
        res += resi
    end
    return res
end

take2(_, v) = v

# calculate res_i = fold(res_i, \sum_{j = 1} A_{i, j} min(u_i, v_j)), i = 1..N
function minplus_bilinear!(foldl::Function,
                           res::AbstractVector{T},
                           A::AbstractMatrix{T},
                           u::AbstractVector{T},
                           v::AbstractVector{T}) where T
    length(res) == length(u) == size(A, 2) ||
        throw(DimensionMismatch("res ($(length(res))), u ($(length(u))) or A $(size(A)) size mismatch"))
    length(v) == size(A, 1) ||
        throw(DimensionMismatch("v ($(length(v))) and A $(size(A)) size mismatch"))
    n = length(res)
    @inbounds for i in axes(A, 2)
        ui = u[i]
        ioffset = (i-1)*size(A, 1)
        x = zero(T)
        if zero(T) < ui < one(T) # if ui == 0, nothing is done
            for j in eachindex(v)
                x += A[ioffset+j] * min(ui, v[j])
            end
        elseif ui == one(T)
            for j in eachindex(v)
                x += A[ioffset+j] * v[j] # FIXME dotprod() faster?
            end
        end
        res[i] = foldl(res[i], x)
    end
    return res
end

function maxplus_linear(v::AbstractVector, A::SparseMatrixCSC)
    length(v) == size(A, 1) ||
        throw(DimensionMismatch("v ($(length(v))) and A $(size(A)) size mismatch"))
    Arows = rowvals(A)
    Avals = nonzeros(A)
    res = 0.0
    @inbounds for i in axes(A, 2)
        wmax = 0.0
        for j in nzrange(A, i)
            jr = Arows[j]
            wmax = max(wmax, v[jr] * Avals[j])
        end
        res += wmax
    end
    return res
end

function maxplus_linear!(res::AbstractVector, A::SparseMatrixCSC, u::AbstractVector)
    length(res) == size(A, 1) ||
        throw(DimensionMismatch("result ($(length(res))) and A $(size(A)) size mismatch"))
    length(u) == size(A, 2) ||
        throw(DimensionMismatch("u ($(length(u))) and A $(size(A)) size mismatch"))
    Arows = rowvals(A)
    Avals = nonzeros(A)
    fill!(res, 0.0)
    @inbounds for i in eachindex(u)
        ((ui = u[i]) > 0) || continue
        for j in nzrange(A, i)
            jr = Arows[j]
            res[jr] = max(res[jr], ui * Avals[j])
        end
    end
    return res
end
