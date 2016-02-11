facts("CoverEnumerator") do # FIXME use weights
    context("[:a] [:b] [::c] [:a :b :c], mask=[:a :b]") do
        sm = OESC.SetMosaic([Set(Symbol[:a]), Set(Symbol[:b]),
                              Set(Symbol[:c]), Set(Symbol[:a, :b, :c])])
        sm_ab = OESC.mask(sm, Set(Symbol[:a, :b]))

        # low prior probability to select sets, high probability to miss active element, so select abc
        cover_etor = OESC.CoverEnumerator(sm_ab, OESC.CoverParams(0.5, 0.1))
        cover_coll = collect(cover_etor; max_set_score=0.0)
        @fact length(cover_coll) --> 1
        cover_coll3 = collect(cover_etor; max_set_score=-1.0)
        @fact length(cover_coll3) --> 0

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then a b c
        cover_etor2 = OESC.CoverEnumerator(sm_ab, OESC.CoverParams(0.1, 0.1))
        cover_coll2 = collect(cover_etor2; max_set_score=0.0)
        @pending length(cover_coll2) --> 2
    end
end
