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

@testset "MaskedSetMosaic" begin
    @testset "empty" begin
        sm = SetMosaic(Set{Symbol}[])
        msm = mask(sm, Set{Symbol}())
        @test unmask(msm) == sm
        @test nelements(msm) == 0
        @test ntiles(msm) == 0
        @test nsets(msm) == 0
        @test nmasked_pertile(msm) == Int[]
        @test nunmasked_pertile(msm) == Int[]
    end

    @testset "empty but with elements" begin
        sm = SetMosaic(Set{Symbol}[], Set([:a, :b]))
        msm = mask(sm, Set([:a]))
        @test unmask(msm) == sm
        @test nelements(msm) == 2
        @test ntiles(msm) == 0
        @test nsets(msm) == 0
        @test nmasked(msm) == 1
        @test nunmasked(msm) == 1
        @test nmasked_pertile(msm) == Int[]
        @test nunmasked_pertile(msm) == Int[]
    end

    @testset "[:a :b] [:c :d] [:a :b :c :d], mask=[:a :b]" begin
        sm = SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c, :d])])

        msm = mask(sm, Set([:a, :b]))

        @test nelements(msm) == 4
        @test nmasked(msm) == 2
        @test nunmasked(msm) == 2
        @test nsets(msm) == 2
        @test length(nmasked_pertile(msm)) == 2
        @test length(nmasked_perset(msm)) == 2
        @test ntiles(msm) == 2
        # FIXME tile indices not stable
        if nmasked_pertile(msm) == Int[2, 0]
            @test_broken tile(msm, 1) == Int[1, 2]
            @test_broken tile(msm, 2) == Int[3, 4]
            @test nmasked_pertile(msm) == Int[2, 0]
            @test nunmasked_pertile(msm) == Int[0, 2]
        else
            @test_broken tile(msm, 1) == Int[3, 4]
            @test_broken tile(msm, 2) == Int[1, 2]
            @test nmasked_pertile(msm) == Int[0, 2]
            @test nunmasked_pertile(msm) == Int[2, 0]
        end
        @test nmasked_perset(msm) == Int[2, 2]
        @test nunmasked_perset(msm) == Int[0, 2]

        msm_copy = copy(msm)
        @test nelements(msm_copy) == nelements(msm)
        @test nmasked(msm_copy) == nmasked(msm)
        @test msm_copy.original == msm.original
        @test nsets(msm_copy) == nsets(msm)

        # mask with nonexisting element
        msm2 = mask(sm, Set([:a, :b, :g]))
        @test nelements(msm2) == 4
        @test ntiles(msm2) == 2
        @test nmasked(msm2) == 2
        @test nunmasked(msm2) == 2
        @test nsets(msm2) == 2

        # mask with max_overlap_logpvalue, [:a :b :c :d] is excluded
        msm3 = mask(sm, Set([:a, :b]), max_overlap_logpvalue=-0.1)
        @test nelements(msm3) == 4
        @test ntiles(msm3) == 2
        @test nmasked(msm3) == 2
        @test nunmasked(msm3) == 2
        @test nsets(msm3) == 1
    end

    @testset "A=[:a :b] B=[:c :d] C=[:a :b :c :d], mask=[:a :b]" begin
        sm = SetMosaic(Dict(:A=>Set([:a, :b]), :B=>Set([:c, :d]), :C=>Set([:a, :b, :c, :d])))
        msm = mask(sm, Set([:a, :b]))

        @test nelements(msm) == 4
        @test nsets(msm) == 2
        @test ntiles(msm) == 2
        # FIXME tile indices not stable
        if nmasked_pertile(msm) == Int[2, 0]
            @test_broken tile(msm, 1) == Int[1, 2]
            @test_broken tile(msm, 2) == Int[3, 4]
            @test nmasked_pertile(msm) == Int[2, 0]
            @test nunmasked_pertile(msm) == Int[0, 2]
        else
            @test_broken tile(msm, 1) == Int[3, 4]
            @test_broken tile(msm, 2) == Int[1, 2]
            @test nmasked_pertile(msm) == Int[0, 2]
            @test nunmasked_pertile(msm) == Int[2, 0]
        end
    end

    @testset "filter!()" begin
        sm = SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c, :d])])
        msm = mask(sm, Set([:b, :c]))

        @test nelements(msm) == 4
        @test ntiles(msm) == 2
        @test nsets(msm) == 3
        @test nmasked_pertile(msm) == Int[1, 1]
        @test nunmasked_pertile(msm) == Int[1, 1]

        filter!(msm, Bool[true, true, false])
        @test nelements(msm) == 4
        @test ntiles(msm) == 2
        @test nsets(msm) == 2
        @test nmasked_pertile(msm) == Int[1, 1]
        @test nunmasked_pertile(msm) == Int[1, 1]

        filter!(msm, Bool[false, true])
        @test nelements(msm) == 4
        @test ntiles(msm) == 1
        @test nsets(msm) == 1
        @test nmasked_pertile(msm) == Int[1]
        @test nunmasked_pertile(msm) == Int[1]
    end
end
