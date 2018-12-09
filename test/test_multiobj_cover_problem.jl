@testset "MultiobjCoverProblem" begin
    @testset "empty" begin
        sm = SetMosaic(Set{Symbol}[])
        msm = mask(sm, [Set{Symbol}()])
        problem = MultiobjCoverProblem(msm)
        @test nvars(problem) == 0
        @test_skip nmasks(problem) == 1
        @test_throws DimensionMismatch rawscore(problem, [1.0])
        @test OESC.miscover_score(problem, Float64[]) == (0.0, 0.0)
        @test rawscore(problem, Float64[]) == (0.0, 0.0, 0.0, 0.0)
        @test score(problem, Float64[]) == (0.0, 0.0)

        res = optimize(problem)
        @test res.weights == Vector{Float64}()
        @test res.agg_total_score == 0.0
    end

    @testset "[a]" begin # FIXME take element detection probability into account
        # empty mask problem
        empty_prob = MultiobjCoverProblem(mask(SetMosaic([Set([:a])]), [Set{Symbol}()]))
        @test nvars(empty_prob) == 0
        @test OESC.miscover_score(empty_prob, Float64[]) == (0.0, 0.0)
        @test rawscore(empty_prob, Float64[]) == (0.0, 0.0, 0.0, 0.0)
        @test score(empty_prob, Float64[]) == (0.0, 0.0)
        empty_res = optimize(empty_prob)
        @test empty_res.weights == Vector{Float64}()
        @test empty_res.agg_total_score == 0.0

        # enabled because sel prob is high, and although :a is all elements it needs to be covered
        prob1_mosaic = mask(SetMosaic([Set([:a])], Set([:a])), [Set([:a])], min_nmasked=1)
        en1_prob = MultiobjCoverProblem(prob1_mosaic, CoverParams(sel_prob=0.5, uncovered_factor=1.0))
        @test nvars(en1_prob) == 1
        @test_throws DimensionMismatch rawscore(en1_prob, Float64[])
        @test OESC.miscover_score(en1_prob, [0.0]) == (1.0, 0.0)
        @test OESC.miscover_score(en1_prob, [1.0]) == (0.0, 0.0)
        en1_rawscore = rawscore(en1_prob, [1.0])
        @test en1_rawscore == (en1_prob.var_scores[1], 0.0, 0.0, 0.0)
        @test rawscore(en1_prob, [0.0]) == (0.0, 0.0, 1.0, 0.0)
        @test en1_rawscore[1] < 1.0 # rawscore[3] at [0.0], required for enabling
        @test rawscore(en1_prob, [0.5]) == (0.5*en1_prob.var_scores[1], 0.0, 0.5, 0.0)
        @test aggscore(en1_prob, [1.0]) < aggscore(en1_prob, [0.9]) < aggscore(en1_prob, [0.5]) < aggscore(en1_prob, [0.0])
        en1_res = optimize(en1_prob, MultiobjOptimizerParams(ϵ=0.01))
        @test en1_res.weights == ones(Float64, 1)

        # disabled because :a is all elements and selection probability low
        dis1_prob = MultiobjCoverProblem(prob1_mosaic, CoverParams(sel_prob=0.01, uncovered_factor=1.0))
        @test nvars(dis1_prob) == 1
        @test dis1_prob.var_scores[1] > en1_prob.var_scores[1]
        dis1_rawscore = rawscore(dis1_prob, [1.0])
        @test dis1_rawscore == (dis1_prob.var_scores[1], 0.0, 0.0, 0.0)
        @test rawscore(dis1_prob, [0.0]) == (0.0, 0.0, 1.0, 0.0)
        @test dis1_rawscore[1] > 1.0 # rawscore[3] at [0.0], required for disabling
        @test aggscore(dis1_prob, [1.0]) > aggscore(dis1_prob, [0.5]) > aggscore(dis1_prob, [0.0])
        dis1_res = optimize(dis1_prob)
        @test dis1_res.weights == zeros(Float64, 1)

        # enabled because there is :a and :b
        en2_prob = MultiobjCoverProblem(mask(SetMosaic([Set([:a])], Set([:a, :b])),
                                             [Set([:a])], min_nmasked=1),
                                        CoverParams(sel_prob=0.9))
        @test nvars(en2_prob) == 1
        en2_res = optimize(en2_prob)
        @test en2_res.weights == ones(Float64, 1)
    end

    @testset "[a b] [c d] [a b c d], mask=[a b]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c])]);
        sm_ab = mask(sm, [Set([:a, :b])], min_nmasked=1)

        problem_def = MultiobjCoverProblem(sm_ab)
        res_def = optimize(problem_def)
        @test res_def.weights ≈ [1.0, 0.0] atol=0.01

        problem_no_penalty = MultiobjCoverProblem(sm_ab, CoverParams(setXset_factor=0.0, covered_factor=0.0, sel_prob=1.0))
        @test nvars(problem_no_penalty) == 2
        res_no_penalty = optimize(problem_no_penalty, MultiobjOptimizerParams(ϵ=[0.01, 0.01]))
        @show res_no_penalty.total_score
        @test res_no_penalty.weights ≈ [1.0, 1.0] atol=0.01
    end

    @testset "[a b] [b c] [a b c], mask=[b]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b]), Set([:b, :c]), Set([:a, :b, :c])])
        sm_b = mask(sm, [Set([:b])], min_nmasked=1)

        prob_hi_sXs = MultiobjCoverProblem(sm_b, CoverParams(setXset_factor=10.0, sel_prob=1E-25))
        @test nvars(prob_hi_sXs) == 3
        res_ignore_overlap = optimize(prob_hi_sXs, MultiobjOptimizerParams(ϵ=0.001))
        @test res_ignore_overlap.weights ≈ [0.0, 0.0, 0.0] atol=0.01

        problem_b = MultiobjCoverProblem(sm_b, CoverParams(setXset_factor=1.0, uncovered_factor=1.0, sel_prob=0.5))
        res_b = optimize(problem_b)
        # [:a :b] is as good as [:b :c]
        @test aggscore(problem_b, [0.0, 1.0, 0.0]) < aggscore(problem_b, [0.0, 0.0, 1.0])
        @test aggscore(problem_b, [0.0, 1.0, 0.0]) < aggscore(problem_b, [0.0, 0.0, 0.0])
        @test aggscore(problem_b, [0.0, 1.0, 0.0]) == aggscore(problem_b, [1.0, 0.0, 0.0])
        @test max(res_b.weights[1], res_b.weights[2]) ≈ 1.0 atol=0.01
        @test res_b.weights[3] ≈ 0.0 atol = 0.01
    end

    @testset "[a b d] [b c d] [c] [d] [a b c d e] [c d e], mask=[a b c]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b, :d]), Set([:b, :c, :d]), Set([:c]), Set([:d]),
                        Set([:a, :b, :c, :d, :e]), Set([:c, :d, :e, :f])],
                        Set([:a, :b, :c, :d, :e, :f]))
        sm_abc = mask(sm, [Set([:a, :b, :c])], min_nmasked=1)

        # low prior probability to select sets, high probability to miss active element, so select abc
        prob_hi_sXs = MultiobjCoverProblem(sm_abc, CoverParams(setXset_factor=10.0, uncovered_factor=1.0, sel_prob=0.1))
        @test nvars(prob_hi_sXs) == 5 # d is out
        @test aggscore(prob_hi_sXs, [0.0, 0.0, 0.0, 1.0, 0.0]) < aggscore(prob_hi_sXs, [0.0, 0.0, 0.0, 0.0, 0.0])
        @test aggscore(prob_hi_sXs, [0.0, 0.0, 0.0, 1.0, 0.0]) < aggscore(prob_hi_sXs, [0.0, 1.0, 0.0, 0.0, 0.0])
        res_ignore_overlap = optimize(prob_hi_sXs, MultiobjOptimizerParams(ϵ=[0.001, 0.001]))
        #@show problem_ab_lowp
        @test aggscore(prob_hi_sXs, [0.0, 0.0, 0.0, 1.0, 0.0]) < aggscore(prob_hi_sXs, [1.0, 0.0, 0.0, 0.0, 0.0])
        @test res_ignore_overlap.weights ≈ [0.0, 0.0, 0.0, 1.0, 0.0]

        # higher prior probability to select sets, lower probability to miss active element, so select a and b
        problem_low_penalty = MultiobjCoverProblem(sm_abc, CoverParams(setXset_factor=1.0, uncovered_factor=1.0, sel_prob=0.9))
        @test aggscore(problem_low_penalty, [1.0, 0.0, 1.0, 0.0, 0.0]) < aggscore(problem_low_penalty, [1.0, 1.0, 0.0, 0.0, 0.0])
        @test aggscore(problem_low_penalty, [1.0, 0.0, 1.0, 0.0, 0.0]) < aggscore(problem_low_penalty, [0.0, 0.0, 0.0, 1.0, 0.0])
        res_low_penalty = optimize(problem_low_penalty, MultiobjOptimizerParams(ϵ=[0.001, 0.001]))
        @test res_low_penalty.weights ≈ [1.0, 0.0, 1.0, 0.0, 0.0]
    end

    @testset "multimask: [a b d] [b c d] [c] [d] [a b c d e] [c d e], mask=[[a b c] [b e]]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b, :d]), Set([:b, :c, :d]), Set([:c]), Set([:d]),
                        Set([:a, :b, :c, :d, :e]), Set([:c, :d, :e, :f])],
                        Set([:a, :b, :c, :d, :e, :f]))
        sm_abc = mask(sm, [Set([:a, :b, :c]), Set([:b, :e])], min_nmasked=1)

        # lower prior probability to select sets, high overlap penalty, so select abd and c FIXME update desc
        prob_hi_sXs = MultiobjCoverProblem(sm_abc, CoverParams(setXset_factor=1.0, uncovered_factor=1.0, covered_factor=0.5, sel_prob=0.6))
        @test_skip nmasks(prob_hi_sXs) == 2
        @test nvars(prob_hi_sXs) == 5 # d is out
        @test aggscore(prob_hi_sXs, [1.0, 0.0, 1.0, 0.0, 0.0]) < aggscore(prob_hi_sXs, [0.0, 0.0, 0.0, 0.0, 0.0])
        @test aggscore(prob_hi_sXs, [0.0, 0.0, 0.0, 1.0, 0.0]) < aggscore(prob_hi_sXs, [0.0, 0.0, 0.0, 0.0, 0.0])
        @test aggscore(prob_hi_sXs, [1.0, 0.0, 1.0, 0.0, 0.0]) > aggscore(prob_hi_sXs, [0.0, 0.0, 0.0, 1.0, 0.0])
        res_hi_sXs = optimize(prob_hi_sXs, MultiobjOptimizerParams(ϵ=[0.001, 0.001]))
        #@show problem_ab_lowp
        @test res_hi_sXs.weights ≈ [0.0, 0.0, 0.0, 1.0, 0.0] atol=1E-2

        # higher prior probability to select sets, no overlap penalty, so select abd, bcde, c and abcde + cdef
        prob_low_sXs = MultiobjCoverProblem(sm_abc, CoverParams(setXset_factor=0.01, covered_factor=1E-5, uncovered_factor=1.0, sel_prob=0.9))
        @test_skip nmasks(prob_low_sXs) == 2
        @test aggscore(prob_low_sXs, [1.0, 0.0, 1.0, 0.0, 0.0]) < aggscore(prob_low_sXs, [0.0, 0.0, 0.0, 0.0, 0.0])
        @test aggscore(prob_low_sXs, [1.0, 0.0, 1.0, 1.0, 0.0]) < aggscore(prob_low_sXs, [1.0, 0.0, 1.0, 0.0, 0.0])
        @test aggscore(prob_low_sXs, [1.0, 0.0, 1.0, 1.0, 0.0]) > aggscore(prob_low_sXs, [0.0, 0.0, 1.0, 1.0, 0.0])
        res_low_sXs = optimize(prob_low_sXs, MultiobjOptimizerParams(ϵ=[0.001, 0.001]))
        #@show res_low_penalty.weights
        @test res_low_sXs.weights ≈ [0.0, 0.0, 1.0, 1.0, 0.0] atol=1E-2# [1.0, 0.0, 1.0, 1.0, 0.0] atol=1E-2
    end
end
