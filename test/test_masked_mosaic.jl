@testset "MaskedSetMosaic" begin
    using OptEnrichedSetCover: MaskOverlap, MaskedSetMosaic

    @testset "empty" begin
        sm = SetMosaic(Set{Symbol}[])
        msm = mask(sm, [Set{Symbol}()])
        @test msm isa MaskedSetMosaic
        @test originalmosaic(msm) == sm
        @test nelements(msm) == 0
        @test nsets(msm) == 0
        @test nmasks(msm) == 1
        @test msm.ix2mask == [1]
    end

    @testset "empty but with elements" begin
        sm = SetMosaic(Set{Symbol}[], Set([:a, :b]))
        msm = mask(sm, [Set([:a])])
        @test originalmosaic(msm) == sm
        @test nelements(msm) == 2
        @test nsets(msm) == 0
        @test nmasks(msm) == 1
        @test nmasked(msm, 1) == 1
        @test nunmasked(msm, 1) == 1
        @test msm.ix2mask == [1]

        # named
        nmsm = mask(sm, Dict(:X => Set([:a])))
        @test originalmosaic(nmsm) == sm
        @test nelements(nmsm) == 2
        @test nsets(nmsm) == 0
        @test nmasks(nmsm) == 1
        @test nmasked(nmsm, :X) == 1
        @test nunmasked(nmsm, :X) == 1
        @test nmsm.ix2mask == [:X]
    end

    @testset "[:a :b] [:c :d] [:a :b :c :d], mask=[:a :b]" begin
        sm = SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c, :d])])

        msm = mask(sm, [Set([:a, :b])])

        @test nelements(msm) == 4
        @test nmasked(msm, 1) == 2
        @test nunmasked(msm, 1) == 2
        @test nsets(msm) == 2
        @test nmasks(msm) == 1
        @test msm.set2masks == Dict(1 => [MaskOverlap(1, 2, 0)],
                                    3 => [MaskOverlap(1, 2, 2)])

        msm_copy = copy(msm)
        @test nelements(msm_copy) == nelements(msm)
        @test nmasks(msm_copy) == nmasks(msm)
        @test msm_copy.original === msm.original
        @test msm_copy.elmasks !== msm.elmasks
        @test msm_copy.elmasks == msm.elmasks
        @test msm_copy.total_masked !== msm.total_masked
        @test msm_copy.total_masked == msm.total_masked
        @test msm_copy.set2masks !== msm.set2masks

        # mask with nonexisting element
        msm2 = mask(sm, [Set([:a, :b, :g])])
        @test nelements(msm2) == 4
        @test nmasked(msm2, 1) == 2
        @test nunmasked(msm2, 1) == 2
        @test nsets(msm2) == 2
        @test msm2.set2masks == Dict(1 => [MaskOverlap(1, 2, 0)],
                                     3 => [MaskOverlap(1, 2, 2)])

        # mask with max_overlap_logpvalue, [:a :b :c :d] is excluded
        msm3 = mask(sm, [Set([:a, :b])], max_overlap_logpvalue=-0.1)
        @test nelements(msm3) == 4
        @test nmasked(msm3, 1) == 2
        @test nunmasked(msm3, 1) == 2
        @test nsets(msm3) == 1
        @test msm3.set2masks == Dict(1 => [MaskOverlap(1, 2, 0)])
    end

    @testset "A=[:a :b] B=[:c :d] C=[:a :b :c :d], mask=[:a :b]" begin
        sm = SetMosaic(Dict(:A=>Set([:a, :b]), :B=>Set([:c, :d]), :C=>Set([:a, :b, :c, :d])))
        msm = mask(sm, [Set([:a, :b])])

        @test nelements(msm) == 4
        @test nsets(msm) == 2
        @test nmasks(msm) == 1
        @test msm.set2masks == Dict(1 => [MaskOverlap(1, 2, 0)],
                                    3 => [MaskOverlap(1, 2, 2)])
    end

    @testset "multimask: A=[:a :b] B=[:c :d] C=[:a :b :c :d], mask=[[:a :b] [:a :b :c] [:a :b :d]]" begin
        sm = SetMosaic(Dict(:A=>Set([:a, :b]), :B=>Set([:c, :d]), :C=>Set([:a, :b, :c, :d])))
        msm = mask(sm, [Set([:a, :b]), Set([:a, :b, :c]), Set([:a, :b, :d])], min_nmasked=1)

        @test nelements(msm) == 4
        @test nsets(msm) == 3
        @test nmasks(msm) == 3
        @test msm.set2masks == Dict(1 => [MaskOverlap(1, 2, 0), MaskOverlap(2, 2, 0), MaskOverlap(3, 2, 0)],
                                    2 => [MaskOverlap(2, 1, 1), MaskOverlap(3, 1, 1)],
                                    3 => [MaskOverlap(1, 2, 2), MaskOverlap(2, 3, 1), MaskOverlap(3, 3, 1)])
        msm2 = mask(sm, [Set([:a, :b]), Set([:a, :b, :c]), Set([:a, :b, :d])], min_nmasked=2)
        @test nsets(msm2) == 2
        @test msm2.set2masks == Dict(1 => [MaskOverlap(1, 2, 0), MaskOverlap(2, 2, 0), MaskOverlap(3, 2, 0)],
                                     3 => [MaskOverlap(1, 2, 2), MaskOverlap(2, 3, 1), MaskOverlap(3, 3, 1)])

        nmsm2 = mask(sm, Dict(:X => Set([:a, :b]), :Y => Set([:a, :b, :c]), :Z => Set([:a, :b, :d])), min_nmasked=2)
        @test nsets(nmsm2) == 2
        @test nmsm2.ix2mask == [:Z, :X, :Y]
        @test nmsm2.set2masks == Dict(1 => [MaskOverlap(1, 2, 0), MaskOverlap(2, 2, 0), MaskOverlap(3, 2, 0)],
                                      3 => [MaskOverlap(1, 3, 1), MaskOverlap(2, 2, 2), MaskOverlap(3, 3, 1)])
    end
end
