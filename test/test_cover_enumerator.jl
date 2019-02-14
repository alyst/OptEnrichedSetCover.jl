global tested_problem_types = [:multiobjective]

@testset "CoverEnumerator problem_type=$problem_type" for problem_type in tested_problem_types # FIXME use weights
    @testset "[a b c abc], mask=[a b]" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])])
        sm_ab = mask(sm, [Set([:a, :b])], min_nmasked=1)

        # low penality to select sets, high probability to miss active element, so select a and b, and abc
        cover_params1 = CoverParams(sel_tax=-log(0.5), uncovered_factor=0.5)
        cover_problem = MultiobjCoverProblem(sm_ab, cover_params1)
        @test nvars(cover_problem) == 3
        @test aggscore([1.0, 1.0, 0.0], cover_problem) < aggscore([0.0, 0.0, 0.0], cover_problem)
        @test aggscore([0.0, 0.0, 1.0], cover_problem) < aggscore([0.0, 0.0, 0.0], cover_problem)
        @test aggscore([1.0, 1.0, 0.0], cover_problem) < aggscore([1.0, 1.0, 1.0], cover_problem)
        @test aggscore([0.0, 0.0, 1.0], cover_problem) < aggscore([1.0, 1.0, 1.0], cover_problem)
        cover_coll1a = collect(sm_ab, cover_params1, CoverEnumerationParams(max_set_score=10.0), problem_type=problem_type)
        @test length(cover_coll1a) == 2
        cover_coll1b = collect(sm_ab, cover_params1, CoverEnumerationParams(max_set_score=0.4), problem_type=problem_type)
        @test length(cover_coll1b) == 1
        cover_coll1c = collect(sm_ab, cover_params1, CoverEnumerationParams(max_set_score=-1.0), problem_type=problem_type)
        @test length(cover_coll1c) == 0

        # higher penalty to select sets, high probability to miss active element, so select abc
        cover_params2 = CoverParams(sel_tax=-log(0.75))
        cover_coll2a = collect(sm_ab, cover_params2, CoverEnumerationParams(max_set_score=10.0), problem_type=problem_type)
        if problem_type==:quadratic # FIXME workfor all
            @test length(cover_coll2a) == 1
        end
        cover_coll2b = collect(sm_ab, cover_params2, CoverEnumerationParams(max_set_score=-1.0), problem_type=problem_type)
        @test length(cover_coll2b) == 0

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then abc
        cover_params3 = CoverParams(sel_tax=0.0)
        cover_coll3 = collect(sm_ab, cover_params3, CoverEnumerationParams(max_set_score=10.0), problem_type=problem_type)
        @test length(cover_coll3) == 2 # FIXME use a and b weights
    end

    @testset "[a b c abc de], mask=[a b]" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                       Set([:a, :b, :c, :d, :e]))

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then abc
        cover_coll = collect(mask(sm, [Set([:a, :b])], min_nmasked=1), CoverParams(sel_tax=0.0),
                             CoverEnumerationParams(max_set_score=10.0), problem_type=problem_type, ϵ=0.01)
        @test length(cover_coll) == 2
        @test best_aggscore(cover_coll.results[1]) <= best_aggscore(cover_coll.results[2])
    end

    @testset "DataFrame(CoverCollection)" begin
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                       Set([:a, :b, :c, :d, :e]))
        sm_ab = mask(sm, [Set([:a, :b])], min_nmasked=1)

        # higher prior probability to select sets, lower probability to miss active element, so select a and b, then a b c
        cover_coll = collect(sm_ab, CoverParams(sel_tax=0.0),
                             CoverEnumerationParams(max_set_score=10.0),
                             problem_type=problem_type, ϵ=0.01)

        df = DataFrame(cover_coll, sm)
        @test size(df, 1) == 3
        @test df.cover_ix == [1, 1, 2]
        @test df.set_id == [1, 2, 4]
        @test df == DataFrame(cover_coll, sm, min_nmasked=1)

        df2 = DataFrame(cover_coll, sm, min_nmasked=0)
        @test size(df2, 1) == 3
        @test df2.cover_ix == [1, 1, 2]
        @test df2.set_id == [1, 2, 4]
    end

    @testset "multimask: [abd bcd c d abcde cde], mask=[abc be]" begin # FIXME take weights into account
        sm = SetMosaic(Dict("abd" => Set([:a, :b, :d]), "bcd" => Set([:b, :c, :d]),
                            "c" => Set([:c]), "d" => Set([:d]),
                            "abcde" => Set([:a, :b, :c, :d, :e]),
                            "cdef" => Set([:c, :d, :e, :f])),
                        Set([:a, :b, :c, :d, :e, :f]))
        sm_abc_be = mask(sm, [Set([:a, :b, :c]), Set([:b, :e])], min_nmasked=1)

        # higher prior probability to select sets, high overlap penalty,
        # high penalty for covering unmasked, so select [abd c], [bcd] and [abcde]
        cover_coll = collect(sm_abc_be, CoverParams(setXset_factor=1.0, covered_factor=0.05, uncovered_factor=0.1, sel_tax=-log(0.9)),
                             CoverEnumerationParams(max_set_score=0.0),
                             problem_type==:quadratic ?
                                OESC.QuadraticOptimizerParams() :
                                OESC.MultiobjOptimizerParams(ϵ=0.01),
                             false)

        df = DataFrame(cover_coll, sm, best_only=true)
        @test size(df, 1) == 8
        @test df.cover_ix == [1, 1, 1, 1, 2, 2, 3, 3]
        @test df.mask_ix == [1, 2, 1, 2, 1, 2, 1, 2]
        @test df.set_id == ["c", "c", "abd", "abd", "abcde", "abcde", "bcd", "bcd"]
    end
end
