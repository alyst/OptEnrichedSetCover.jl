@testset "CoverEnumerator" begin # FIXME use weights
    @testset "[:a] [:b] [:c] [:a :b :c], mask=[:a :b]" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])])
        sm_ab = mask(sm, Set([:a, :b]))

        # low prior probability to select sets, high probability to miss active element, so select abc
        cover_etor = CoverEnumerator(sm_ab, CoverParams(0.5, 0.1))
        cover_coll = collect(cover_etor; max_set_score=0.0)
        @test_broken length(cover_coll) == 1
        cover_coll3 = collect(cover_etor; max_set_score=-1.0)
        @test length(cover_coll3) == 0

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then a b c
        cover_etor2 = CoverEnumerator(sm_ab, CoverParams(0.1, 0.1))
        cover_coll2 = collect(cover_etor2; max_set_score=0.0, max_cover_score_delta=0.0)
        @test_broken length(cover_coll2) == 2
    end

    @testset "[:a] [:b] [:c] [:a :b :c] :d :e, mask=[:a :b]" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                            Set([:a, :b, :c, :d, :e]))
        sm_ab = mask(sm, Set(Symbol[:a, :b]))

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then a b c
        cover_etor2 = CoverEnumerator(sm_ab, CoverParams(0.1, 0.1))
        cover_coll2 = collect(cover_etor2; max_set_score=0.0, max_cover_score_delta=0.0)
        @test_broken length(cover_coll2) == 2
        @test_broken cover_coll2.variants[1].score <= cover_coll2.variants[2].score
    end

    @testset "convert(DataFrame)" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                            Set([:a, :b, :c, :d, :e]))
        sm_ab = mask(sm, Set(Symbol[:a, :b]))

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then a b c
        cover_etor2 = CoverEnumerator(sm_ab, CoverParams(0.1, 0.1))
        cover_coll2 = collect(cover_etor2; max_set_score=0.0, max_cover_score_delta=0.0)

        df = convert(DataFrame, cover_coll2, sm)
        @test_broken size(df, 1) == 3
        @test_broken df[:cover_ix] == [1, 1, 2]
        @test_broken df[:set_id] == [1, 2, 4]
    end
end
