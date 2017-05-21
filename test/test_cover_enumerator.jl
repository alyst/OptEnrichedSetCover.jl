@testset "CoverEnumerator" begin # FIXME use weights
    @testset "[:a] [:b] [:c] [:a :b :c], mask=[:a :b]" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])])
        sm_ab = mask(sm, Set([:a, :b]))

        # low penality to select sets, high probability to miss active element, so select abc
        cover_params1 = CoverParams(sel_prob=1.0)
        cover_coll1a = collect(sm_ab, cover_params1, CoverEnumerationParams(max_set_score=0.0))
        @test length(cover_coll1a) == 1
        cover_coll1b = collect(sm_ab, cover_params1, CoverEnumerationParams(max_set_score=-1.0))
        @test length(cover_coll1b) == 0

        # high penalty to select sets, high probability to miss active element, so select abc
        cover_params2 = CoverParams(sel_prob=0.01)
        cover_coll2a = collect(sm_ab, cover_params2, CoverEnumerationParams(max_set_score=0.0))
        @test length(cover_coll2a) == 1
        cover_coll2b = collect(sm_ab, cover_params2, CoverEnumerationParams(max_set_score=-1.0))
        @test length(cover_coll2b) == 0

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then abc
        cover_params3 = CoverParams(sel_prob=1.0)
        cover_coll3 = collect(sm_ab, cover_params3, CoverEnumerationParams(max_set_score=0.0, max_cover_score_delta=0.0))
        #@show cover_coll2
        @test_broken length(cover_coll3) == 2 # FIXME use a and b weights
    end

    @testset "[:a] [:b] [:c] [:a :b :c] :d :e, mask=[:a :b]" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                            Set([:a, :b, :c, :d, :e]))

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then abc
        cover_coll = collect(mask(sm, Set(Symbol[:a, :b])), CoverParams(sel_prob=1.0),
                             CoverEnumerationParams(max_set_score=0.0, max_cover_score_delta=0.0))
        @test length(cover_coll) == 2
        @test cover_coll.variants[1].score <= cover_coll.variants[2].score
    end

    @testset "DataFrame(CoverCollection)" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                            Set([:a, :b, :c, :d, :e]))
        sm_ab = mask(sm, Set(Symbol[:a, :b]))

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then a b c
        cover_coll = collect(sm_ab, CoverParams(sel_prob=1.0),
                             CoverEnumerationParams(max_set_score=0.0, max_cover_score_delta=0.0))

        df = DataFrame(cover_coll, sm)
        @test size(df, 1) == 3
        # FIXME mess with DataFrames Nullables, should switch to DataTables at some point
        @test_skip df[:cover_ix] == [1, 1, 2]
        @test_skip df[:set_id] == [1, 2, 4]
    end
end
