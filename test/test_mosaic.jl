@testset "SetMosaic" begin
    elms(sm::SetMosaic, ixs::AbstractVector) =
        getindex.(Ref(sm.ix2elm), ixs)
    tile_elms(sm::SetMosaic, ix::Integer) =
        elms(sm, tile(sm, ix))

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

    @testset "[] [], no elements" begin
        sm = SetMosaic([Set{Symbol}(), Set{Symbol}()])

        @test nelements(sm) == 0
        @test ntiles(sm) == 0
        @test nsets(sm) == 1
        @test sm.set2ix == Dict(1=>1, 2=>1)
        @test sm.ix2set == [1]
    end

    @testset "[:a]" begin
        sm = SetMosaic([Set([:a])])

        @test nelements(sm) == 1
        @test ntiles(sm) == 1
        @test tile_elms(sm, 1) == [:a]
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
        @test tile_elms(sm, 1) == [:b]
        @test tile_elms(sm, 2) == [:a]
        @test nsets(sm) ==2
    end

    @testset "[:a] [:a]" begin
        sm = SetMosaic([Set([:a]), Set([:a])])

        @test nelements(sm) == 1
        @test ntiles(sm) == 1
        @test tile_elms(sm, 1) == [:a]
        @test nsets(sm) == 1
        @test sm.set2ix == Dict(1=>1, 2=>1)
        @test sm.ix2set == [1]
    end

    @testset "[:a] [:b] [:a :b]" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:a, :b])])

        @test nelements(sm) == 2
        @test ntiles(sm) == 2
        @test tile_elms(sm, 1) == [:a]
        @test tile_elms(sm, 2) == [:b]
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
        @test tile_elms(sm, 1) == [:a, :b]
        @test tile_elms(sm, 2) == [:d]
        @test tile_elms(sm, 3) == [:c]
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

    @testset "[:a] [:a] [:a :b] [:b :c] [:b c] [] []" begin
        sm = SetMosaic([Set([:a]), Set([:a]), Set([:a, :b]), Set([:b, :c]), Set([:b, :c]), Set{Symbol}()],
                       Set([:a, :b, :c, :d]),
                       [0.5, 1.0, 0.7, 0.3, 0.1, 0.2])

        @test nelements(sm) == 4
        @test ntiles(sm) == 4
        @test nsets(sm) == 4
        @test sm.set2ix == Dict(1=>1, 2=>1, 3=>2, 4=>3, 5=>3, 6=>4)
        @test sm.ix2set == [2, 3, 4, 6]
        @test sm.set_relevances == [1.0, 0.7, 0.3, 0.2]
    end
end
