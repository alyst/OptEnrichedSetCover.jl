facts("SparseMaskMatrix") do
    context("empty") do
        mask = OESC.SparseMaskMatrix()
        @fact size(mask) --> (0, 0)
        @fact size(mask, 1) --> 0
        @fact size(mask, 2) --> 0
        @fact_throws size(mask, 0) ArgumentError
        @fact_throws size(mask, 3) ArgumentError

        empty2 = OESC.SparseMaskMatrix(2, 3)
        @fact size(empty2) --> (2, 3)
        @fact slice(empty2, :, 1) --> Int[]
        @fact slice(empty2, :, 2) --> Int[]
        @fact slice(empty2, :, 3) --> Int[]
    end

    context("1x1") do
        empty = OESC.SparseMaskMatrix(1, 1, [1, 1], Int[])
        @fact size(empty) --> (1, 1)
        @fact size(empty, 1) --> 1
        @fact size(empty, 2) --> 1
        @fact empty[:, 1] --> Int[]
        @fact slice(empty, :, 1) --> Int[]

        nonempty = OESC.SparseMaskMatrix(1, 1, [1, 2], Int[1])
        @fact size(nonempty) --> (1, 1)
        @fact size(nonempty, 1) --> 1
        @fact size(nonempty, 2) --> 1
        @fact nonempty[:, 1] --> [1]
        @fact slice(nonempty, :, 1) --> [1]
    end

    context("sparse_mask(Collection{Set}, elm2ix)") do
        sm = OESC.sparse_mask([Set([:a, :b]), Set{Symbol}(), Set([:c, :d]), Set([:a, :c])],
                         Dict(:a=>1, :b=>2, :c=>3, :d=>4, :e=>5))
        @fact size(sm) --> (5, 4)
        @fact slice(sm, :, 1) --> [1, 2]
        @fact slice(sm, :, 2) --> Int[]
        @fact slice(sm, :, 3) --> [3, 4]
        @fact slice(sm, :, 4) --> [1, 3]
    end

    context("sparse_mask(Vector{Vector{Int}})") do
        sm = OESC.sparse_mask(4, 3, Vector{Int}[[1,2], Int[], [3,4]])
        @fact size(sm) --> (4, 3)
        @fact slice(sm, :, 1) --> [1, 2]
        @fact slice(sm, :, 2) --> Int[]
        @fact slice(sm, :, 3) --> [3, 4]
    end

    context("convert(Matrix{Bool})") do
        sm = convert(OESC.SparseMaskMatrix, [false false true; true false false])
        @fact size(sm) --> (2, 3)
        @fact slice(sm, :, 1) --> [2]
        @fact slice(sm, :, 2) --> Int[]
        @fact slice(sm, :, 3) --> [1]
    end
end
