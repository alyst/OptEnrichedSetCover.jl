@testset "MaskedSetMosaic" begin
    using OptEnrichedSetCover: MaskedSet, MaskedSetMosaic

    @testset "MaskedSet" begin
        @test setsize(MaskedSet(1, 1, 2, 0)) == 2
        @test setsize(MaskedSet(1, 1, 3, 2)) == 5
    end

    @testset "empty" begin
        sm = SetMosaic(Set{Symbol}[])
        msm = mask(sm, [Set{Symbol}()])
        @test msm isa MaskedSetMosaic
        @test unmask(msm) == sm
        @test nelements(msm) == 0
        @test nsets(msm) == 0
        @test nmaskedsets(msm) == 0
        @test nmasks(msm) == 1
    end

    @testset "empty but with elements" begin
        sm = SetMosaic(Set{Symbol}[], Set([:a, :b]))
        msm = mask(sm, [Set([:a])])
        @test unmask(msm) == sm
        @test nelements(msm) == 2
        @test nsets(msm) == 0
        @test nmaskedsets(msm) == 0
        @test nmasks(msm) == 1
        @test nmasked(msm, 1) == 1
        @test nunmasked(msm, 1) == 1
    end

    @testset "[:a :b] [:c :d] [:a :b :c :d], mask=[:a :b]" begin
        sm = SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c, :d])])

        msm = mask(sm, [Set([:a, :b])])

        @test nelements(msm) == 4
        @test nmasked(msm, 1) == 2
        @test nunmasked(msm, 1) == 2
        @test nsets(msm) == 2
        @test nmaskedsets(msm) == 2
        @test nmasks(msm) == 1
        @test maskedset(msm, 1) == MaskedSet(1, 1, 2, 0)
        @test maskedset(msm, 2) == MaskedSet(1, 3, 2, 2)

        msm_copy = copy(msm)
        @test nelements(msm_copy) == nelements(msm)
        @test nmasks(msm_copy) == nmasks(msm)
        @test msm_copy.original === msm.original
        @test msm_copy.elmasks !== msm.elmasks
        @test msm_copy.elmasks == msm.elmasks
        @test msm_copy.total_masked !== msm.total_masked
        @test msm_copy.total_masked == msm.total_masked
        @test msm_copy.maskedsets !== msm.maskedsets
        @test msm_copy.maskedsets == msm.maskedsets

        # mask with nonexisting element
        msm2 = mask(sm, [Set([:a, :b, :g])])
        @test nelements(msm2) == 4
        @test nmasked(msm2, 1) == 2
        @test nunmasked(msm2, 1) == 2
        @test nsets(msm2) == 2
        @test nmaskedsets(msm2) == 2
        @test maskedset(msm, 1) == MaskedSet(1, 1, 2, 0)
        @test maskedset(msm, 2) == MaskedSet(1, 3, 2, 2)

        # mask with max_overlap_logpvalue, [:a :b :c :d] is excluded
        msm3 = mask(sm, [Set([:a, :b])], max_overlap_logpvalue=-0.1)
        @test nelements(msm3) == 4
        @test nmasked(msm3, 1) == 2
        @test nunmasked(msm3, 1) == 2
        @test nsets(msm3) == 1
        @test nmaskedsets(msm3) == 1
        @test maskedset(msm, 1) == MaskedSet(1, 1, 2, 0)
    end

    @testset "A=[:a :b] B=[:c :d] C=[:a :b :c :d], mask=[:a :b]" begin
        sm = SetMosaic(Dict(:A=>Set([:a, :b]), :B=>Set([:c, :d]), :C=>Set([:a, :b, :c, :d])))
        msm = mask(sm, [Set([:a, :b])])

        @test nelements(msm) == 4
        @test nsets(msm) == 2
        @test nmaskedsets(msm) == 2
        @test nmasks(msm) == 1
        @test maskedset(msm, 1) == MaskedSet(1, 1, 2, 0)
        @test maskedset(msm, 2) == MaskedSet(1, 3, 2, 2)
    end

    @testset "multimask: A=[:a :b] B=[:c :d] C=[:a :b :c :d], mask=[[:a :b] [:a :b :c] [:a :b :d]]" begin
        sm = SetMosaic(Dict(:A=>Set([:a, :b]), :B=>Set([:c, :d]), :C=>Set([:a, :b, :c, :d])))
        msm = mask(sm, [Set([:a, :b]), Set([:a, :b, :c]), Set([:a, :b, :d])], min_nmasked=1)

        @test nelements(msm) == 4
        @test nmaskedsets(msm) == 8
        @test nsets(msm) == 3
        @test nmasks(msm) == 3
        @test maskedset(msm, 1) == MaskedSet(1, 1, 2, 0)
        @test maskedset(msm, 2) == MaskedSet(1, 3, 2, 2)
        @test maskedset(msm, 3) == MaskedSet(2, 1, 2, 0)
        @test maskedset(msm, 4) == MaskedSet(2, 2, 1, 1)
        @test maskedset(msm, 5) == MaskedSet(2, 3, 3, 1)
        @test maskedset(msm, 6) == MaskedSet(3, 1, 2, 0)
        @test maskedset(msm, 7) == MaskedSet(3, 2, 1, 1)
        @test maskedset(msm, 8) == MaskedSet(3, 3, 3, 1)

        msm2 = mask(sm, [Set([:a, :b]), Set([:a, :b, :c]), Set([:a, :b, :d])], min_nmasked=2)
        @test nmaskedsets(msm2) == 6
        @test nsets(msm2) == 2
    end

    @testset "filter!()" begin
        sm = SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c, :d])])
        msm = mask(sm, [Set([:b, :c])], min_nmasked=1)

        @test nelements(msm) == 4
        @test nmaskedsets(msm) == 3
        @test nsets(msm) == 3
        @test maskedset(msm, 1) == MaskedSet(1, 1, 1, 1)
        @test maskedset(msm, 2) == MaskedSet(1, 2, 1, 1)
        @test maskedset(msm, 3) == MaskedSet(1, 3, 2, 2)

        filter!(msm, Bool[true, true, false])
        @test nelements(msm) == 4
        @test nmaskedsets(msm) == 2
        @test_broken nsets(msm) == 2
        @test maskedset(msm, 1) == MaskedSet(1, 1, 1, 1)
        @test maskedset(msm, 2) == MaskedSet(1, 2, 1, 1)

        filter!(msm, Bool[false, true])
        @test nelements(msm) == 4
        @test nmaskedsets(msm) == 1
        @test_broken nsets(msm) == 1
        @test maskedset(msm, 1) == MaskedSet(1, 2, 1, 1)
    end
end
