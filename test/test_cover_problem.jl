@testset "CoverProblem" begin
    @testset "CoverParams" begin
        @test_throws ArgumentError CoverParams(sel_prob=-0.1)
        @test_throws ArgumentError CoverParams(sel_prob=1.1)
    end

    @testset "empty" begin
        sm = SetMosaic(Set{Symbol}[])
        msm = mask(sm, [Set{Symbol}()])
        problem = CoverProblem(msm)
        @test nvars(problem) == 0
        @test_skip nmasks(problem) == 1
        @test score(problem, Float64[]) == 0.0

        res = optimize(problem)
        @test res.weights == Vector{Float64}()
        @test res.total_score == 0.0
    end

    @testset "[a]" begin # FIXME take element detection probability into account
        # empty mask problem
        empty_problem = CoverProblem(mask(SetMosaic([Set([:a])]), [Set{Symbol}()]))
        @test nvars(empty_problem) == 0
        empty_res = optimize(empty_problem)
        @test empty_res.weights == Vector{Float64}(0)
        @test empty_res.total_score == 0.0

        # disabled because :a is all elements
        dis_problem = CoverProblem(mask(SetMosaic([Set([:a])], Set([:a])), [Set{Symbol}([:a])], min_nmasked=1))
        @test nvars(dis_problem) == 1
        dis_res = optimize(dis_problem)
        @test dis_res.weights == zeros(Float64, 1)

        # enabled because there is :a and :b
        en_problem = CoverProblem(mask(SetMosaic([Set([:a])], Set([:a, :b])), [Set{Symbol}([:a])], min_nmasked=1))
        @test nvars(en_problem) == 1
        en_res = optimize(en_problem)
        @test en_res.weights == ones(Float64, 1)
    end

    @testset "[a b] [c d] [a b c d], mask=[a b]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c])]);
        sm_ab = mask(sm, [Set(Symbol[:a, :b])])

        problem_def = CoverProblem(sm_ab)
        res_def = optimize(problem_def)
        @test res_def.weights == [1.0, 0.0]

        problem_no_penalty = CoverProblem(sm_ab, CoverParams(setXset_factor=0.0, sel_prob=1.0))
        @test nvars(problem_no_penalty) == 2
        res_no_penalty = optimize(problem_no_penalty)
        @test res_no_penalty.weights == [1.0, 1.0]
    end

    @testset "[a b] [b c] [a b c], mask=[b]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b]), Set([:b, :c]), Set([:a, :b, :c])])
        sm_b = mask(sm, [Set([:b])], min_nmasked=1)

        problem_ignore_overlap = CoverProblem(sm_b, CoverParams(setXset_factor=10.0, sel_prob=1E-25))
        @test nvars(problem_ignore_overlap) == 3
        res_ignore_overlap = optimize(problem_ignore_overlap)
        @test res_ignore_overlap.weights ≈ [0.0, 0.0, 0.0]# atol=1E-4

        problem_b = CoverProblem(sm_b)
        res_b = optimize(problem_b)
        # [:a :b] is as good as [:b :c]
        @test res_b.weights[1] ≈ res_b.weights[2]
        @test res_b.weights[3] ≈ 0.0
    end

    @testset "[a b d] [b c d] [c] [d] [a b c d e] [c d e], mask=[a b c]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b, :d]), Set([:b, :c, :d]), Set([:c]), Set([:d]),
                        Set([:a, :b, :c, :d, :e]), Set([:c, :d, :e, :f])],
                        Set([:a, :b, :c, :d, :e, :f]))
        sm_abc = mask(sm, [Set([:a, :b, :c])], min_nmasked=1)

        # low prior probability to select sets, high probability to miss active element, so select abc
        problem_ignore_overlap = CoverProblem(sm_abc, CoverParams(setXset_factor=10.0, sel_prob=0.1))
        @test nvars(problem_ignore_overlap) == 5 # d is out
        @test score(problem_ignore_overlap, [0.0, 0.0, 0.0, 1.0, 0.0]) < score(problem_ignore_overlap, [1.0, 1.0, 0.0, 0.0, 0.0])
        res_ignore_overlap = optimize(problem_ignore_overlap)
        #@show problem_ab_lowp
        @test_broken res_ignore_overlap.weights ≈ [0.0, 0.0, 0.0, 1.0, 0.0]# atol=1E-4

        # higher prior probability to select sets, lower probability to miss active element, so select a and b
        problem_low_penalty = CoverProblem(sm_abc, CoverParams(setXset_factor=1.0, sel_prob=0.9))
        @test score(problem_low_penalty, [1.0, 0.0, 1.0, 0.0, 0.0]) < score(problem_low_penalty, [1.0, 1.0, 0.0, 0.0, 0.0])
        @test score(problem_low_penalty, [1.0, 0.0, 1.0, 0.0, 0.0]) < score(problem_low_penalty, [0.0, 0.0, 0.0, 1.0, 0.0])
        res_low_penalty = optimize(problem_low_penalty)
        @test find(res_low_penalty.weights) == [1, 3]
    end

    @testset "multimask: [a b d] [b c d] [c] [d] [a b c d e] [c d e], mask=[[a b c] [b e]]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b, :d]), Set([:b, :c, :d]), Set([:c]), Set([:d]),
                        Set([:a, :b, :c, :d, :e]), Set([:c, :d, :e, :f])],
                        Set([:a, :b, :c, :d, :e, :f]))
        sm_abc = mask(sm, [Set([:a, :b, :c]), Set([:b, :e])], min_nmasked=1)

        # lower prior probability to select sets, high overlap penalty, so select abd and c
        problem_ignore_overlap = CoverProblem(sm_abc, CoverParams(setXset_factor=1.0, sel_prob=0.6))
        @test_skip nmasks(problem_ignore_overlap) == 2
        @test nvars(problem_ignore_overlap) == 9 # c and 2xd are out
        res_ignore_overlap = optimize(problem_ignore_overlap)
        #@show problem_ab_lowp
        @test res_ignore_overlap.weights ≈ [1.0, 0.0, 1.0, #=d,=# 0.0, 0.0,
                                            0.0, 0.0, #=c,=# #=d,=# 0.0, 0.0] atol=1E-4

        # higher prior probability to select sets, no overlap penalty, so select abd, bcde, c and abcde + cdef
        problem_low_penalty = CoverProblem(sm_abc, CoverParams(setXset_factor=0.05, sel_prob=0.9))
        res_low_penalty = optimize(problem_low_penalty)
        #@show res_low_penalty.weights
        @test find(res_low_penalty.weights) == [1, 2, 3, 4, 8]
    end
end
