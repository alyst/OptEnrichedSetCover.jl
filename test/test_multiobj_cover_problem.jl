@testset "MultiobjCoverProblem" begin
    @testset "empty" begin
        sm = SetMosaic(Set{Symbol}[])
        msm = mask(sm, [Set{Symbol}()])
        problem = MultiobjCoverProblem(msm)
        @test nvars(problem) == 0
        @test_skip nmasks(problem) == 1
        @test_throws DimensionMismatch score([1.0], problem)
        @test OESC.miscover_score(Float64[], problem) == (0.0, 0.0)
        @test score(Float64[], problem) == (0.0, 0.0, 0.0, 0.0)
        @test aggscore(Float64[], problem) == 0.0

        res = optimize(problem)
        @test isempty(res)
        @test nvars(res) == 0
        @test nsolutions(res) == 0
        @test best_index(res) == 0
    end

    @testset "[a]" begin # FIXME take element detection probability into account
        # empty mask problem
        empty_prob = MultiobjCoverProblem(mask(SetMosaic([Set([:a])]), [Set{Symbol}()]))
        @test nvars(empty_prob) == 0
        @test OESC.miscover_score(Float64[], empty_prob) == (0.0, 0.0)
        @test score(Float64[], empty_prob) == (0.0, 0.0, 0.0, 0.0)
        @test aggscore(Float64[], empty_prob) == 0.0
        empty_res = optimize(empty_prob)
        @test isempty(empty_res)
        @test nvars(empty_res) == 0
        @test nsolutions(empty_res) == 0

        # enabled because sel prob is high, and although :a is all elements it needs to be covered
        prob1_mosaic = mask(SetMosaic([Set([:a])], Set([:a])), [Set([:a])], min_nmasked=1)
        en1_prob = MultiobjCoverProblem(prob1_mosaic, CoverParams(sel_tax=-log(0.5), uncovered_factor=1.0))
        @test nvars(en1_prob) == 1
        @test_throws DimensionMismatch score(Float64[], en1_prob)
        @test OESC.miscover_score([0.0], en1_prob) == (1.0, 0.0)
        @test OESC.miscover_score([1.0], en1_prob) == (0.0, 0.0)
        en1_score = score([1.0], en1_prob)
        @test en1_score == (en1_prob.var_scores[1], 0.0, 0.0, 0.0)
        @test score([0.0], en1_prob) == (0.0, 0.0, 1.0, 0.0)
        @test en1_score[1] < 1.0 # score[3] at [0.0], required for enabling
        @test score([0.5], en1_prob) == (0.5*en1_prob.var_scores[1], 0.0, 0.5, 0.0)
        @test aggscore([1.0], en1_prob) < aggscore([0.9], en1_prob) <
              aggscore([0.5], en1_prob) < aggscore([0.0], en1_prob)
        en1_res = optimize(en1_prob, MultiobjOptimizerParams(ϵ=0.001))
        @test nsolutions(en1_res) > 0
        @test nvars(en1_res) == 1
        @test best_varweights(en1_res) == ones(Float64, 1)

        # disabled because :a is all elements and selection probability low
        dis1_prob = MultiobjCoverProblem(prob1_mosaic, CoverParams(sel_tax=-log(0.01), uncovered_factor=1.0))
        @test nvars(dis1_prob) == 1
        @test dis1_prob.var_scores[1] > en1_prob.var_scores[1]
        dis1_score = score([1.0], dis1_prob)
        @test dis1_score == (dis1_prob.var_scores[1], 0.0, 0.0, 0.0)
        @test score([0.0], dis1_prob) == (0.0, 0.0, 1.0, 0.0)
        @test dis1_score[1] > 1.0 # score[3] at [0.0], required for disabling
        @test aggscore([1.0], dis1_prob) > aggscore([0.5], dis1_prob) > aggscore([0.0], dis1_prob)
        dis1_res = optimize(dis1_prob)
        @test nsolutions(dis1_res) > 0
        @test nvars(dis1_res) == 1
        @test best_varweights(dis1_res) == zeros(Float64, 1)

        # enabled because there is :a and :b
        en2_prob = MultiobjCoverProblem(mask(SetMosaic([Set([:a])], Set([:a, :b])),
                                             [Set([:a])], min_nmasked=1),
                                        CoverParams(sel_tax=-log(0.9)))
        @test nvars(en2_prob) == 1
        en2_res = optimize(en2_prob)
        @test nsolutions(en2_res) > 0
        @test nvars(en2_res) == 1
        @test best_varweights(en2_res) == ones(Float64, 1)
    end

    @testset "[a b] [c d] [a b c d], mask=[a b]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c])])
        sm_ab = mask(sm, [Set([:a, :b])], min_nmasked=1)

        problem_def = MultiobjCoverProblem(sm_ab) # cd set is dropped

        @test OESC.var_scores([0.0, 0.0], problem_def) == (0.0, 0.0)
        @test OESC.var_scores([1.0, 0.0], problem_def)[2] == 0.0
        @test OESC.var_scores([0.0, 1.0], problem_def)[2] == 0.0
        @test OESC.var_scores([1.0, 1.0], problem_def)[2] > 0

        @test OESC.miscover_score([0.0, 0.0], problem_def) == (2.0, 0.0)
        @test OESC.miscover_score([0.25, 0.0], problem_def) == (1.5, 0.0)
        @test OESC.miscover_score([0.75, 0.0], problem_def) == (0.5, 0.0)
        @test OESC.miscover_score([1.0, 0.0], problem_def) == (0.0, 0.0)
        @test OESC.miscover_score([0.0, 1.0], problem_def) == (0.0, 1.0)
        @test OESC.miscover_score([1.0, 1.0], problem_def) == (0.0, 1.0)

        res_def = optimize(problem_def)
        @test nsolutions(res_def) > 0
        @test best_varweights(res_def) ≈ [1.0, 0.0] atol=0.01

        problem_no_penalty = MultiobjCoverProblem(sm_ab, CoverParams(setXset_factor=0.0, covered_factor=0.0, sel_tax=0.0))
        @test nvars(problem_no_penalty) == 2
        res_no_penalty = optimize(problem_no_penalty, MultiobjOptimizerParams(ϵ=[0.001, 0.001], WeightDigits=nothing))
        @test nvars(res_no_penalty) == 2
        @test best_varweights(res_no_penalty) ≈ [1.0, 1.0] atol=0.01
    end

    @testset "[a b] [b c] [a b c], mask=[b]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b]), Set([:b, :c]), Set([:a, :b, :c])])
        sm_b = mask(sm, [Set([:b])], min_nmasked=1)

        prob_hi_sXs = MultiobjCoverProblem(sm_b, CoverParams(setXset_factor=10.0, sel_tax=-log(1E-25)))
        @test nvars(prob_hi_sXs) == 3
        res_ignore_overlap = optimize(prob_hi_sXs, MultiobjOptimizerParams(ϵ=0.001))
        @test nvars(res_ignore_overlap) == 3
        @test best_varweights(res_ignore_overlap) ≈ [0.0, 0.0, 0.0] atol=0.01

        problem_b = MultiobjCoverProblem(sm_b, CoverParams(setXset_factor=1.0, uncovered_factor=1.0, sel_tax=-log(0.5)))
        res_b = optimize(problem_b)
        @test nsolutions(res_b) > 0
        @test nvars(res_b) == 3
        # [:a :b] is as good as [:b :c]
        @test aggscore([0.0, 1.0, 0.0], problem_b) < aggscore([0.0, 0.0, 1.0], problem_b)
        @test aggscore([0.0, 1.0, 0.0], problem_b) < aggscore([0.0, 0.0, 0.0], problem_b)
        @test aggscore([0.0, 1.0, 0.0], problem_b) == aggscore([1.0, 0.0, 0.0], problem_b)
        @test max(best_varweights(res_b)[1], best_varweights(res_b)[2]) ≈ 1.0 atol=0.01
        @test best_varweights(res_b)[3] ≈ 0.0 atol = 0.01

        @test OESC.miscover_score([1.0, 0.0, 0.0], problem_b) == (0.0, 1.0)
        @test OESC.miscover_score([0.0, 1.0, 0.0], problem_b) == (0.0, 1.0)
        @test OESC.miscover_score([0.0, 0.5, 0.0], problem_b) == (0.5, 0.5)
        @test OESC.miscover_score([0.0, 0.0, 1.0], problem_b) == (0.0, 2.0)
        @test OESC.miscover_score([1.0, 0.0, 1.0], problem_b) == (0.0, 2.0)
        @test OESC.miscover_score([1.0, 1.0, 1.0], problem_b) == (0.0, 2.0)
        @test OESC.miscover_score([1.0, 1.0, 0.5], problem_b) == (0.0, 2.0)
    end

    @testset "[a b d] [b c d] [c] [d] [a b c d e] [c d e], mask=[a b c]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b, :d]), Set([:b, :c, :d]), Set([:c]), Set([:d]),
                        Set([:a, :b, :c, :d, :e]), Set([:c, :d, :e, :f])],
                        Set([:a, :b, :c, :d, :e, :f]))
        sm_abc = mask(sm, [Set([:a, :b, :c])], min_nmasked=1)

        # low prior probability to select sets, high probability to miss active element, so select abc
        prob_hi_sXs = MultiobjCoverProblem(sm_abc, CoverParams(setXset_factor=10.0, uncovered_factor=1.0, sel_tax=-log(0.1)))
        @test nvars(prob_hi_sXs) == 5 # d is out

        @test aggscore([0.0, 0.0, 0.0, 1.0, 0.0], prob_hi_sXs) < aggscore([0.0, 0.0, 0.0, 0.0, 0.0], prob_hi_sXs)
        @test aggscore([0.0, 0.0, 0.0, 1.0, 0.0], prob_hi_sXs) < aggscore([0.0, 1.0, 0.0, 0.0, 0.0], prob_hi_sXs)
        res_ignore_overlap = optimize(prob_hi_sXs, MultiobjOptimizerParams(ϵ=[0.001, 0.001]))
        @test aggscore([0.0, 0.0, 0.0, 1.0, 0.0], prob_hi_sXs) < aggscore([1.0, 0.0, 0.0, 0.0, 0.0], prob_hi_sXs)
        @test best_varweights(res_ignore_overlap) ≈ [0.0, 0.0, 0.0, 1.0, 0.0]

        # higher prior probability to select sets, lower probability to miss active element, so select a and b
        problem_low_penalty = MultiobjCoverProblem(sm_abc, CoverParams(setXset_factor=1.0, uncovered_factor=1.0, sel_tax=-log(0.9)))
        @test aggscore([1.0, 0.0, 1.0, 0.0, 0.0], problem_low_penalty) < aggscore([1.0, 1.0, 0.0, 0.0, 0.0], problem_low_penalty)
        @test aggscore([1.0, 0.0, 1.0, 0.0, 0.0], problem_low_penalty) < aggscore([0.0, 0.0, 0.0, 1.0, 0.0], problem_low_penalty)
        res_low_penalty = optimize(problem_low_penalty, MultiobjOptimizerParams(ϵ=[0.001, 0.001]))
        @test best_varweights(res_low_penalty) ≈ [1.0, 0.0, 1.0, 0.0, 0.0]
    end

    @testset "multimask: [a b d] [b c d] [c] [d] [a b c d e] [c d e], mask=[[a b c] [b e]]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b, :d]), Set([:b, :c, :d]), Set([:c]), Set([:d]),
                        Set([:a, :b, :c, :d, :e]), Set([:c, :d, :e, :f])],
                        Set([:a, :b, :c, :d, :e, :f]))
        sm_abc = mask(sm, [Set([:a, :b, :c]), Set([:b, :e])], min_nmasked=1)

        # lower prior probability to select sets, high overlap penalty, so select abd and c FIXME update desc
        prob_hi_sXs = MultiobjCoverProblem(sm_abc, CoverParams(setXset_factor=1.0, uncovered_factor=1.0, covered_factor=0.5, sel_tax=-log(0.6)))
        @test_skip nmasks(prob_hi_sXs) == 2
        @test nvars(prob_hi_sXs) == 5 # d is out: abd bcd c abcde cdef

        @test OESC.miscover_score([0.0, 0.0, 0.0, 0.0, 0.0], prob_hi_sXs) == (5.0, 0.0)
        @test OESC.miscover_score([1.0, 0.0, 0.0, 0.0, 0.0], prob_hi_sXs) == (2.0, 1.0 + 0.9 * 2.0)
        @test OESC.miscover_score([0.0, 1.0, 0.0, 0.0, 0.0], prob_hi_sXs) == (2.0, 1.0 + 0.9 * 2.0)
        @test OESC.miscover_score([0.0, 0.0, 1.0, 0.0, 0.0], prob_hi_sXs) == (4.0, 0.0 + 0.9 * 1.0)

        @test aggscore([1.0, 0.0, 1.0, 0.0, 0.0], prob_hi_sXs) < aggscore([0.0, 0.0, 0.0, 0.0, 0.0], prob_hi_sXs)
        @test aggscore([0.0, 0.0, 0.0, 1.0, 0.0], prob_hi_sXs) < aggscore([0.0, 0.0, 0.0, 0.0, 0.0], prob_hi_sXs)
        @test aggscore([1.0, 0.0, 1.0, 0.0, 0.0], prob_hi_sXs) > aggscore([0.0, 0.0, 0.0, 1.0, 0.0], prob_hi_sXs)
        res_hi_sXs = optimize(prob_hi_sXs, MultiobjOptimizerParams(ϵ=[0.001, 0.001]))
        @test best_varweights(res_hi_sXs) ≈ [0.0, 0.0, 0.0, 1.0, 0.0] atol=1E-2

        # higher prior probability to select sets, no overlap penalty, so select abd, bcd, c and abcde + cdef
        prob_low_sXs = MultiobjCoverProblem(sm_abc, CoverParams(setXset_factor=0.01, covered_factor=1E-5, uncovered_factor=1.0, sel_tax=-log(0.9)))
        @test_skip nmasks(prob_low_sXs) == 2
        @test aggscore([1.0, 0.0, 1.0, 0.0, 0.0], prob_low_sXs) < aggscore([0.0, 0.0, 0.0, 0.0, 0.0], prob_low_sXs)
        @test aggscore([1.0, 0.0, 1.0, 1.0, 0.0], prob_low_sXs) < aggscore([1.0, 0.0, 1.0, 0.0, 0.0], prob_low_sXs)
        @test aggscore([1.0, 0.0, 1.0, 1.0, 0.0], prob_low_sXs) < aggscore([0.0, 0.0, 1.0, 1.0, 0.0], prob_low_sXs)
        @test aggscore([1.0, 1.0, 1.0, 1.0, 0.0], prob_low_sXs) < aggscore([0.0, 0.0, 1.0, 1.0, 0.0], prob_low_sXs)
        @test aggscore([1.0, 1.0, 1.0, 1.0, 0.0], prob_low_sXs) < aggscore([1.0, 0.0, 1.0, 1.0, 0.0], prob_low_sXs)
        res_low_sXs = optimize(prob_low_sXs, MultiobjOptimizerParams(ϵ=[0.001, 0.001]))
        @test best_varweights(res_low_sXs) ≈ [1.0, 1.0, 1.0, 1.0, 0.0] atol=1E-2# [1.0, 0.0, 1.0, 1.0, 0.0] atol=1E-2
    end

    @testset "WeightedSetMosaic-based problem" begin
        @testset "empty" begin
            sm = SetMosaic(Set{Symbol}[])
            wsm = assignweights(sm, [Dict{Int, Float64}()])
            problem = MultiobjCoverProblem(wsm)
            @test nvars(problem) == 0
            @test_skip nexperiments(problem) == 1
            @test_throws DimensionMismatch score([1.0], problem)
            @test score(Float64[], problem) == (0.0, 0.0, 0.0, 0.0)
            @test aggscore(Float64[], problem) == 0.0

            res = optimize(problem)
            @test isempty(res)
            @test nvars(res) == 0
            @test nsolutions(res) == 0
            @test best_index(res) == 0
        end

        @testset "[a]" begin # FIXME take element detection probability into account
            # empty mask problem
            empty_prob = MultiobjCoverProblem(assignweights(SetMosaic([Set([:a])]), [Dict{Int, Float64}()]))
            @test nvars(empty_prob) == 0
            @test score(Float64[], empty_prob) == (0.0, 0.0, 0.0, 0.0)
            @test aggscore(Float64[], empty_prob) == 0.0
            empty_res = optimize(empty_prob)
            @test isempty(empty_res)
            @test nvars(empty_res) == 0
            @test nsolutions(empty_res) == 0

            # enabled because sel prob is high, and although :a is all elements it needs to be covered
            prob1_mosaic = assignweights(SetMosaic([Set([:a])], Set([:a])), [Dict(1 => -1.0)])
            en1_prob = MultiobjCoverProblem(prob1_mosaic, CoverParams(sel_tax=-log(0.5), uncovered_factor=1.0))
            @test nvars(en1_prob) == 1
            @test_throws DimensionMismatch score(Float64[], en1_prob)
            en1_score = score([1.0], en1_prob)
            @test en1_score == (en1_prob.var_scores[1], 0.0, 0.0, 0.0)
            @test score([0.0], en1_prob) == (0.0, 0.0, 1.0, 0.0)
            @test en1_score[1] < 1.0 # score[3] at [0.0], required for enabling
            @test score([0.5], en1_prob) == (0.5*en1_prob.var_scores[1], 0.0, 0.5, 0.0)
            @test aggscore([1.0], en1_prob) < aggscore([0.9], en1_prob) <
                aggscore([0.5], en1_prob) < aggscore([0.0], en1_prob)
            en1_res = optimize(en1_prob, MultiobjOptimizerParams(ϵ=0.001))
            @test nsolutions(en1_res) > 0
            @test nvars(en1_res) == 1
            @test best_varweights(en1_res) == ones(Float64, 1)

            # disabled because the weight is zero
            dis1_prob = MultiobjCoverProblem(prob1_mosaic, CoverParams(sel_tax=-log(0.01), uncovered_factor=1.0))
            @test nvars(dis1_prob) == 1
            @test dis1_prob.var_scores[1] > en1_prob.var_scores[1]
            dis1_score = score([1.0], dis1_prob)
            @test dis1_score == (dis1_prob.var_scores[1], 0.0, 0.0, 0.0)
            @test score([0.0], dis1_prob) == (0.0, 0.0, 1.0, 0.0)
            @test dis1_score[1] > 1.0 # score[3] at [0.0], required for disabling
            @test aggscore([1.0], dis1_prob) > aggscore([0.5], dis1_prob) > aggscore([0.0], dis1_prob)
            dis1_res = optimize(dis1_prob)
            @test nsolutions(dis1_res) > 0
            @test nvars(dis1_res) == 1
            @test best_varweights(dis1_res) == zeros(Float64, 1)

            # enabled because the weight is negative
            prob2_mosaic = assignweights(SetMosaic([Set([:a])], Set([:a])), [Dict(1 =>  -1)])
            en2_prob = MultiobjCoverProblem(prob2_mosaic, CoverParams(sel_tax=-log(0.9)))
            @test nvars(en2_prob) == 1
            en2_res = optimize(en2_prob)
            @test nsolutions(en2_res) > 0
            @test nvars(en2_res) == 1
            @test best_varweights(en2_res) == ones(Float64, 1)
        end

        @testset "[a b] [c d] [a b c d], mask=[a b]" begin # FIXME take weights into account
            sm = SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c])])
            sm_ab = assignweights(sm, [Dict(1 => -5, 3 => -3)])

            problem_def = MultiobjCoverProblem(sm_ab) # cd set is dropped

            @test OESC.var_scores([0.0, 0.0], problem_def) == (0.0, 0.0)
            @test OESC.var_scores([1.0, 0.0], problem_def)[2] == 0.0
            @test OESC.var_scores([0.0, 1.0], problem_def)[2] == 0.0
            @test OESC.var_scores([1.0, 1.0], problem_def)[2] > 0

            res_def = optimize(problem_def)
            @test nsolutions(res_def) > 0
            @test best_varweights(res_def) ≈ [1.0, 0.0] atol=0.01

            problem_no_penalty = MultiobjCoverProblem(sm_ab, CoverParams(setXset_factor=0.0, covered_factor=0.0, sel_tax=0.0))
            @test nvars(problem_no_penalty) == 2
            res_no_penalty = optimize(problem_no_penalty, MultiobjOptimizerParams(ϵ=[0.001, 0.001], WeightDigits=nothing))
            @test nvars(res_no_penalty) == 2
            @test best_varweights(res_no_penalty) ≈ [1.0, 1.0] atol=0.01
        end
    end
end
