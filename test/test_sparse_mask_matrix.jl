@testset "SparseMaskMatrix" begin
    using OptEnrichedSetCover.SparseMaskMatrix
    @testset "empty" begin
        mask = SparseMaskMatrix()
        @test size(mask) == (0, 0)
        @test size(mask, 1) == 0
        @test size(mask, 2) == 0
        @test_throws ArgumentError size(mask, 0)
        @test_throws ArgumentError size(mask, 3)

        empty2 = SparseMaskMatrix(2, 3)
        @test size(empty2) == (2, 3)
        @test view(empty2, :, 1) == Int[]
        @test view(empty2, :, 2) == Int[]
        @test view(empty2, :, 3) == Int[]
        @test mask != empty2
        @test empty2 == SparseMaskMatrix(2, 3)
    end

    @testset "1x1" begin
        @inferred SparseMaskMatrix(1, 1, [1, 1], Int[])
        empty = SparseMaskMatrix(1, 1, [1, 1], Int[])
        @test size(empty) == (1, 1)
        @test size(empty, 1) == 1
        @test size(empty, 2) == 1
        @test empty[:, 1] == Int[]
        @test view(empty, :, 1) == Int[]

        nonempty = SparseMaskMatrix(1, 1, [1, 2], Int[1])
        @test size(nonempty) == (1, 1)
        @test size(nonempty, 1) == 1
        @test size(nonempty, 2) == 1
        @test nonempty[:, 1] == [1]
        @test view(nonempty, :, 1) == [1]
        @test nonempty != empty
    end

    @testset "SparseMaskMatrix(Collection{Set}, elm2ix)" begin
        sm = SparseMaskMatrix([Set([:a, :b]), Set{Symbol}(), Set([:c, :d]), Set([:a, :c])],
                              Dict(:a=>1, :b=>2, :c=>3, :d=>4, :e=>5))
        @test size(sm) == (5, 4)
        @test view(sm, :, 1) == [1, 2]
        @test view(sm, :, 2) == Int[]
        @test view(sm, :, 3) == [3, 4]
        @test view(sm, :, 4) == [1, 3]
    end

    @testset "SparseMaskMatrix(Vector{Vector{Int}})" begin
        sm = SparseMaskMatrix(4, Vector{Int}[[1,2], Int[], [3,4]])
        @test size(sm) == (4, 3)
        @test view(sm, :, 1) == [1, 2]
        @test view(sm, :, 2) == Int[]
        @test view(sm, :, 3) == [3, 4]
    end

    @testset "SparseMaskMatrix(Matrix{Bool})" begin
        sm = SparseMaskMatrix([false false true; true false false])
        @test size(sm) == (2, 3)
        @test view(sm, :, 1) == [2]
        @test view(sm, :, 2) == Int[]
        @test view(sm, :, 3) == [1]
        @test sm == SparseMaskMatrix(2, [[2], Int[], [1]])
        @test sm != SparseMaskMatrix(3, [[2], Int[], [1]])
        @test sm != SparseMaskMatrix(2, [[1], Int[], [1]])
        @test sm != SparseMaskMatrix(2, [[2], [1], Int[]])
    end

    @testset "convert(Matrix{Bool}, SparseMaskMatrix)" begin
        sm1 = SparseMaskMatrix()
        @inferred convert(Matrix{Bool}, sm1)
        mtx1 = convert(Matrix{Bool}, sm1)
        @test mtx1 isa Matrix{Bool}
        @test mtx1 == Matrix{Bool}(0, 0)

        sm2 = SparseMaskMatrix([false false true; true false false])
        mtx2 = convert(Matrix{Bool}, sm2)
        @test mtx2 == [false false true; true false false]

        mtx3 = convert(Matrix, sm2)
        @test mtx3 isa Matrix{Bool}
        @test mtx3 == mtx2
    end

    @testset "getindex(SparseMaskMatrix)" begin
        sm = SparseMaskMatrix([false false true; true false false])

        @testset "SMS[:, ?]" begin
            @inferred sm[:, 1]
            @test sm[:, 1] == [2]
            @test sm[:, 2] == Int[]

            @inferred sm[:, [1, 2]]
            sm1 = sm[:, [1, 2]]
            @test sm1 isa SparseMaskMatrix
            @test convert(Matrix, sm1) == [false false; true false]
            @test size(sm1) == (2, 2)
            @inferred sm[:, [false, true, true]]
            @test sm[:, [2, 3]] == sm[:, [false, true, true]]
            @test sm[:, [2, 1, 3, 1]] == SparseMaskMatrix([false false true false; false true false true])
        end

        @testset "SMS[?, :]" begin
            @test_broken sm[1, :] == [false false true] # not implemented yet
            @test sm[[1], :] == SparseMaskMatrix([false false true])
            @inferred sm[[2,1,2], :]
            @test sm[[2,1,2], :] == SparseMaskMatrix([true false false; false false true; true false false])
        end

        @testset "SMS[?, ?]" begin
            @inferred sm[1, 2]
            @test sm[1, 1] == false
            @test sm[1, 3] == true
            @inferred sm[[2,1,2], [1]]
            @test sm[[2,1,2], [1]] == SparseMaskMatrix(reshape([true, false, true], 3, 1))
        end
    end

    @testset "permutedims(SMS, ?)" begin
        @test transpose(SparseMaskMatrix()) == SparseMaskMatrix()
        @test transpose(SparseMaskMatrix(2, 3)) == SparseMaskMatrix(3, 2)

        sm = SparseMaskMatrix([false false true; true false false])
        @inferred permutedims(sm, [1, 2])
        @test permutedims(sm, [1, 2]) == sm
        @inferred permutedims(sm, [2, 1])
        tsm = permutedims(sm, [2, 1])
        @test size(tsm) == (3, 2)
        @test tsm == SparseMaskMatrix([false true; false false; true false])
        @test transpose(sm) == tsm
        @test transpose(transpose(sm)) == sm
    end
end
