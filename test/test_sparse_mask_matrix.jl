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
    end

    @testset "1x1" begin
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
    end
end
