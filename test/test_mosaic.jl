@testset "SetMosaic" begin
    @testset "no sets, no elements" begin
        sm = SetMosaic(Set{Symbol}[])

        @test nelements(sm) == 0
        @test ntiles(sm) == 0
        @test nsets(sm) == 0
    end

    @testset "no sets, :a :b" begin
        sm = SetMosaic(Set{Symbol}[], Set([:a, :b]))
        @test nelements(sm) == 2
        @test ntiles(sm) == 0
        @test nsets(sm) == 0
        @test_throws BoundsError setsize(sm, 0)
        @test_throws BoundsError setsize(sm, 1)
    end

    @testset "[], no elements" begin
        sm = SetMosaic([Set{Symbol}()])

        @test nelements(sm) == 0
        @test ntiles(sm) == 0
        @test nsets(sm) == 1
    end

    @testset "[:a]" begin
        sm = SetMosaic([Set([:a])])

        @test nelements(sm) == 1
        @test ntiles(sm) == 1
        @test tile(sm, 1) == [1]
        @test nsets(sm) == 1
        @test setsize(sm, 1) == 1
    end

    @testset "[:a] []" begin
        sm = SetMosaic([Set([:a]), Set{Symbol}()])

        @test nelements(sm) == 1
        @test ntiles(sm) == 1
        @test tile(sm, 1) == [1]
        @test nsets(sm) == 2
        @test setsize(sm, 1) == 1
        @test setsize(sm, 2) == 0
    end

    @testset "[:a] [:b]" begin
        sm = SetMosaic([Set([:a]), Set([:b])])

        @test nelements(sm) == 2
        @test ntiles(sm) == 2
        @test_broken tile(sm, 1) == [1]
        @test_broken tile(sm, 2) == [2]
        @test nsets(sm) ==2
    end

    @testset "[:a] [:a]" begin
        sm = SetMosaic([Set([:a]), Set([:a])])

        @test nelements(sm) == 1
        @test ntiles(sm) == 1
        @test tile(sm, 1) == [1]
        @test nsets(sm) == 2
    end

    @testset "[:a] [:b] [:a :b]" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:a, :b])])

        @test nelements(sm) == 2
        @test ntiles(sm) == 2
        @test tile(sm, 1) == [1]
        @test tile(sm, 2) == [2]
        @test setsize(sm, 1) == 1
        @test setsize(sm, 2) == 1
        @test setsize(sm, 3) == 2
        @test_throws BoundsError setsize(sm, 4)
        @test nsets(sm) == 3
    end

    @testset "[:a :b] [:c :d] [:a :b :c]" begin
        sm = SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c])])

        @test nelements(sm) == 4
        @test ntiles(sm) == 3
        @test_broken tile(sm, 1) == [1, 2]
        @test_broken tile(sm, 2) == [3]
        @test_broken tile(sm, 3) == [4]
        @test nsets(sm) == 3
    end

    @testset "[:a :b] [:c :d] [:a :b :c :d]" begin
        sm = SetMosaic([Set([:a, :b]), Set([:c, :d]),
                        Set([:a, :b, :c, :d])])

        @test nelements(sm) == 4
        @test ntiles(sm) == 2
        @test tile(sm, 1) == [1, 2]
        @test tile(sm, 2) == [3, 4]
        @test nsets(sm) == 3
    end
end
