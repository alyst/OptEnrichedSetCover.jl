@testset "WeightedSetMosaic" begin
    using OptEnrichedSetCover: WeightedSetMosaic

    @testset "empty" begin
        sm = SetMosaic(Set{Symbol}[])
        wsm = assignweights(sm, [Dict{Int, Float64}()])
        @test wsm isa WeightedSetMosaic
        @test originalmosaic(wsm) == sm
        @test nelements(wsm) == 0
        @test nsets(wsm) == 0
        @test nexperiments(wsm) == 1
    end

    @testset "empty but with elements" begin
        sm = SetMosaic(Set{Symbol}[], Set([:a, :b]))
        wsm = assignweights(sm, [Dict{Int, Float64}()])
        @test wsm isa WeightedSetMosaic
        @test originalmosaic(wsm) == sm
        @test nelements(wsm) == 2
        @test nsets(wsm) == 0
        @test nexperiments(wsm) == 1

        # named
        nwsm = assignweights(sm, Dict(:X => Dict{Int, Float64}()))
        @test originalmosaic(nwsm) == sm
        @test nelements(nwsm) == 2
        @test nsets(nwsm) == 0
        @test nexperiments(nwsm) == 1
    end

    @testset "[:a :b] [:c :d] [:a :b :c :d]" begin
        sm = SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c, :d])])

        wsm = assignweights(sm, [Dict(1 => -5.0, 3 => -2.0)])

        @test nelements(wsm) == 4
        @test nexperiments(wsm) == 1
        @test nsets(wsm) == 2
        @test setweight(wsm, 1, 1) == -5
        @test setweight(wsm, 2, 1) === missing
        @test setweight(wsm, 3, 1) == -2
        @test_throws KeyError setweight(wsm, 1, 2)
        @test_throws KeyError setweight(wsm, 4, 1)

        wsm_copy = copy(wsm)
        @test nelements(wsm_copy) == nelements(wsm)
        @test nexperiments(wsm_copy) == nexperiments(wsm)
        @test wsm_copy.original === wsm.original
        @test wsm_copy.setXexp_weights !== wsm.setXexp_weights
        @test wsm_copy.setXexp_weights == wsm.setXexp_weights

        # ignore weight of nonexisting set
        wsm2 = assignweights(sm, [Dict(4 => -8.0)])
        @test nelements(wsm2) == 4
        @test nsets(wsm2) == 0
        @test_throws KeyError setweight(wsm2, 4, 1)

        # mask with maweightsx_min_weight, [:a :b :c :d] is excluded
        wsm3 = assignweights(sm, [Dict(1 => -5.0, 3 => -2.0)], max_min_weight=-2.5)
        @test nelements(wsm3) == 4
        @test nsets(wsm3) == 1
        @test setweight(wsm3, 1, 1) == -5
        @test setweight(wsm3, 2, 1) === missing
        @test setweight(wsm3, 3, 1) === missing

        wsm4 = assignweights(sm, [Dict(1 => -5.0, 3 => -2.0)], max_setsize=3)
        @test nelements(wsm4) == 4
        @test nsets(wsm4) == 1
        @test setweight(wsm4, 1, 1) == -5
        @test setweight(wsm4, 2, 1) === missing
        @test setweight(wsm4, 3, 1) === missing
    end

    @testset "A=[:a :b] B=[:c :d] C=[:a :b :c :d]" begin
        sm = SetMosaic(Dict(:A=>Set([:a, :b]), :B=>Set([:c, :d]), :C=>Set([:a, :b, :c, :d])))
        wsm = assignweights(sm, [Dict(:A => -5.0, :C => -2.0)])

        @test nelements(wsm) == 4
        @test nsets(wsm) == 2
        @test nexperiments(wsm) == 1
        @test setweight(wsm, :A, 1) == -5
        @test setweight(wsm, :B, 1) === missing
        @test setweight(wsm, :C, 1) == -2
    end

    @testset "multimask: A=[:a :b] B=[:c :d] C=[:a :b :c :d], mask=[[:a :b] [:a :b :c] [:a :b :d]]" begin
        sm = SetMosaic(Dict(:A=>Set([:a, :b]), :B=>Set([:c, :d]), :C=>Set([:a, :b, :c, :d])))
        setweights = [Dict(:A => -5.0,             :C => -2.0),
                      Dict(:A => -6.0, :B => -1.0, :C => -8.0),
                      Dict(            :B => -2.0, :C => -3.0)]
        wsm = assignweights(sm, setweights)
        @test nelements(wsm) == 4
        @test nsets(wsm) == 3
        @test nexperiments(wsm) == 3
        @test_throws KeyError setweight(wsm, :A, 4)
        @test_throws KeyError setweight(wsm, :D, 1)
        @test setweight(wsm, :A, 1) == -5
        @test setweight(wsm, :A, 2) == -6
        @test setweight(wsm, :A, 3) === missing
        @test setweight(wsm, :B, 1) === missing
        @test setweight(wsm, :B, 2) == -1
        @test setweight(wsm, :B, 3) == -2
        @test setweight(wsm, :C, 1) == -2
        @test setweight(wsm, :C, 2) == -8
        @test setweight(wsm, :C, 3) == -3

        wsm2 = assignweights(sm, setweights, max_min_weight=-3)
        @test nsets(wsm2) == 2
        @test setweight(wsm2, :A, 1) == -5
        @test setweight(wsm2, :A, 2) == -6
        @test setweight(wsm2, :A, 3) === missing
        @test setweight(wsm2, :B, 1) === missing
        @test setweight(wsm2, :B, 2) === missing
        @test setweight(wsm2, :B, 3) === missing
        @test setweight(wsm2, :C, 1) == -2
        @test setweight(wsm2, :C, 2) == -8
        @test setweight(wsm2, :C, 3) == -3

        nsetweights = Dict(:X => setweights[1],
                           :Y => setweights[2],
                           :Z => setweights[3])
        nwsm2 = assignweights(sm, nsetweights, max_weight=-2, max_setsize=3)
        @test nsets(nwsm2) == 2 # C filtered out because it's too big
        @test nexperiments(nwsm2) == 3
        @test_throws KeyError setweight(nwsm2, :A, :U)
        @test_throws KeyError setweight(nwsm2, :D, :X)
        @test setweight(nwsm2, :A, :X) == -5
        @test setweight(nwsm2, :A, :Y) == -6
        @test setweight(nwsm2, :A, :Z) === missing
        @test setweight(nwsm2, :B, :X) === missing
        @test setweight(nwsm2, :B, :Y) === missing
        @test setweight(nwsm2, :B, :Z) == -2
        @test setweight(nwsm2, :C, :X) === missing
        @test setweight(nwsm2, :C, :Y) === missing
        @test setweight(nwsm2, :C, :Z) === missing
    end
end
