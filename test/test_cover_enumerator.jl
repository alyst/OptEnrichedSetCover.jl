facts("CoverEnumerator") do # FIXME use weights
    context("[:a] [:b] [:c] [:a :b :c], mask=[:a :b]") do
        sm = OESC.SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])])
        sm_ab = OESC.mask(sm, Set([:a, :b]))

        # low prior probability to select sets, high probability to miss active element, so select abc
        cover_etor = OESC.CoverEnumerator(sm_ab, OESC.CoverParams(0.5, 0.1))
        cover_coll = collect(cover_etor; max_set_score=0.0)
        @fact length(cover_coll) --> 1
        cover_coll3 = collect(cover_etor; max_set_score=-1.0)
        @fact length(cover_coll3) --> 0

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then a b c
        cover_etor2 = OESC.CoverEnumerator(sm_ab, OESC.CoverParams(0.1, 0.1))
        cover_coll2 = collect(cover_etor2; max_set_score=0.0, max_cover_score_delta=0.0)
        @pending length(cover_coll2) --> 2
    end

    context("[:a] [:b] [:c] [:a :b :c] :d :e, mask=[:a :b]") do
        sm = OESC.SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                            Set([:a, :b, :c, :d, :e]))
        sm_ab = OESC.mask(sm, Set(Symbol[:a, :b]))

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then a b c
        cover_etor2 = OESC.CoverEnumerator(sm_ab, OESC.CoverParams(0.1, 0.1))
        cover_coll2 = collect(cover_etor2; max_set_score=0.0, max_cover_score_delta=0.0)
        @fact length(cover_coll2) --> 2
        @fact cover_coll2.variants[1].score --> less_than_or_equal(cover_coll2.variants[2].score)
    end

    context("convert(DataFrame)") do
        sm = OESC.SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                            Set([:a, :b, :c, :d, :e]))
        sm_ab = OESC.mask(sm, Set(Symbol[:a, :b]))

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then a b c
        cover_etor2 = OESC.CoverEnumerator(sm_ab, OESC.CoverParams(0.1, 0.1))
        cover_coll2 = collect(cover_etor2; max_set_score=0.0, max_cover_score_delta=0.0)

        df = convert(DataFrame, cover_coll2, sm)
        @fact size(df, 1) --> 3
        @fact df[:cover_ix] --> [1, 1, 2]
        @fact df[:set_id] --> [1, 2, 4]
    end
end
