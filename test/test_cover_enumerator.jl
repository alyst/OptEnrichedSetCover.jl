@testset "CoverEnumerator" begin # FIXME use weights
    @testset "[:a] [:b] [:c] [:a :b :c], mask=[:a :b]" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])])
        sm_ab = mask(sm, [Set([:a, :b])])

        # low penality to select sets, high probability to miss active element, so select abc
        cover_params1 = CoverParams(sel_prob=1.0)
        cover_coll1a = collect(sm_ab, cover_params1, CoverEnumerationParams(max_set_score=10.0))
        @test length(cover_coll1a) == 1
        cover_coll1b = collect(sm_ab, cover_params1, CoverEnumerationParams(max_set_score=-1.0))
        @test length(cover_coll1b) == 0

        # higher penalty to select sets, high probability to miss active element, so select abc
        cover_params2 = CoverParams(sel_prob=0.75)
        cover_coll2a = collect(sm_ab, cover_params2, CoverEnumerationParams(max_set_score=10.0))
        @test length(cover_coll2a) == 1
        cover_coll2b = collect(sm_ab, cover_params2, CoverEnumerationParams(max_set_score=-1.0))
        @test length(cover_coll2b) == 0

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then abc
        cover_params3 = CoverParams(sel_prob=1.0)
        cover_coll3 = collect(sm_ab, cover_params3, CoverEnumerationParams(max_set_score=10.0))
        #@show cover_coll2
        @test_broken length(cover_coll3) == 2 # FIXME use a and b weights
    end

    @testset "[:a] [:b] [:c] [:a :b :c] :d :e, mask=[:a :b]" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                            Set([:a, :b, :c, :d, :e]))

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then abc
        cover_coll = collect(mask(sm, [Set(Symbol[:a, :b])]), CoverParams(sel_prob=1.0),
                             CoverEnumerationParams(max_set_score=10.0))
        @test length(cover_coll) == 2
        @test cover_coll.results[1].total_score <= cover_coll.results[2].total_score
    end

    @testset "DataFrame(CoverCollection)" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                            Set([:a, :b, :c, :d, :e]))
        sm_ab = mask(sm, [Set(Symbol[:a, :b])])

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then a b c
        cover_coll = collect(sm_ab, CoverParams(sel_prob=1.0),
                             CoverEnumerationParams(max_set_score=10.0))

        df = DataFrame(cover_coll, sm)
        @test size(df, 1) == 3
        @test df[:cover_ix] == [1, 1, 2]
        @test df[:set_id] == [1, 2, 4]
        @test df == DataFrame(cover_coll, sm, report=:covered)

        df2 = DataFrame(cover_coll, sm, report=:matrix)
        @test size(df2, 1) == 3
        @test df2[:cover_ix] == [1, 1, 2]
        @test df2[:set_id] == [1, 2, 4]

        @test_throws ArgumentError DataFrame(cover_coll, sm, report=:unknown)
    end

    @testset "multimask: [a b d] [b c d] [c] [d] [a b c d e] [c d e], mask=[[a b c] [b e]]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b, :d]), Set([:b, :c, :d]), Set([:c]), Set([:d]),
                        Set([:a, :b, :c, :d, :e]), Set([:c, :d, :e, :f])],
                        Set([:a, :b, :c, :d, :e, :f]))
        sm_abc_be = mask(sm, [Set([:a, :b, :c]), Set([:b, :e])])

        # higher prior probability to select sets, no overlap penalty, so select abd, bcde, c and abcde + cdef
        cover_coll = collect(sm_abc_be, CoverParams(setXset_factor=0.05, sel_prob=0.9),
                             CoverEnumerationParams(max_set_score=0.0))

        df = DataFrame(cover_coll, sm, report=:covered)
        @test size(df, 1) == 5
        @test df[:cover_ix] == [1, 1, 1, 1, 1]
        @test df[:mask_ix] == [1, 1, 1, 1, 2]
        @test df[:set_id] == [1, 2, 3, 5, 5]

        df2 = DataFrame(cover_coll, sm, report=:matrix)
        @test size(df2, 1) == 8
        @test df2[:cover_ix] == [1, 1, 1, 1, 0, 0, 0, 1]
        @test df2[:mask_ix] == [1, 1, 1, 1, 2, 2, 2, 2]
        @test df2[:set_id] == [1, 2, 3, 5, 1, 2, 3, 5]
    end
end
