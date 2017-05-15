@testset "CoverProblem" begin
    @testset "CoverParams" begin
    end

    @testset "empty" begin
        sm = SetMosaic(Set{Symbol}[])
        msm = mask(sm, Set{Symbol}())
        problem = CoverProblem(msm)
        @test nsets(problem) == 0
        @test score(problem, Float64[]) == 0.0

        res = optimize(problem)
        @test res.weights == Float64[]
        @test res.score == 0.0
    end

    @testset "[:a]" begin # FIXME take element detection probability into account
        # empty mask problem
        empty_problem = CoverProblem(mask(SetMosaic([Set([:a])]), Set{Symbol}()))
        @test nsets(empty_problem) == 0
        empty_res = optimize(empty_problem)
        @test empty_res.weights == Float64[]
        @test empty_res.score == 0.0

        # disabled because :a is all elements
        dis_problem = CoverProblem(mask(SetMosaic([Set([:a])], Set([:a])), Set{Symbol}([:a])))
        @test nsets(dis_problem) == 1
        dis_res = optimize(dis_problem)
        @test dis_res.weights == [0.0]

        # enabled because there is :a and :b
        en_problem = CoverProblem(mask(SetMosaic([Set([:a])], Set([:a, :b])), Set{Symbol}([:a])))
        @test nsets(en_problem) == 1
        en_res = optimize(en_problem)
        @test en_res.weights == [1.0]
    end

    @testset "[:a :b] [:c :d] [:a :b :c :d], mask=[:a :b]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c, :d])]);
        sm_ab = mask(sm, Set(Symbol[:a, :b]))

        problem_ab = CoverProblem(sm_ab)
        res_ab = optimize(problem_ab)
        @test res_ab.weights ≈ [1.0, 0.0] atol=1E-4

        # FIXME when weights would be available
        problem_ab_lowp = CoverProblem(sm_ab, CoverParams(a=0.1, b=0.1))
        @test nsets(problem_ab_lowp) == 2
        res_ab_lowp = optimize(problem_ab_lowp)
        @test_skip res_ab_lowp.weights ≈ [1.0, 0.0] atol=1E-4
    end

    @testset "[:a :b] [:b :c] [:a :b :c], mask=[:b]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a, :b]), Set([:b, :c]), Set([:a, :b, :c])])
        sm_b = mask(sm, Set([:b]))

        problem_b_lowp = CoverProblem(sm_b, CoverParams(a=0.1, b=0.1))
        @test nsets(problem_b_lowp) == 3
        res_b_lowp = optimize(problem_b_lowp)
        @test_skip res_b_lowp.weights ≈ [0.0, 0.0, 0.0] atol=1E-4

        problem_b = CoverProblem(sm_b)
        res_b = optimize(problem_b)
        # [:a :b] is as good as [:b :c]
        @test res_b.weights[1] ≈ res_b.weights[2]
        @test res_b.weights[3] ≈ 0.0
    end

    @testset "[:a] [:b] [::c] [:a :b :c], mask=[:a :b]" begin # FIXME take weights into account
        sm = SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                            Set([:a, :b, :c, :d, :e]))
        sm_ab = mask(sm, Set([:a, :b]))

        # low prior probability to select sets, high probability to miss active element, so select abc
        problem_ab_lowp = CoverProblem(sm_ab, CoverParams(a=0.4, b=0.1))
        @test nsets(problem_ab_lowp) == 3 # c is out
        @test_skip score(problem_ab_lowp, [0.0, 0.0, 1.0]) < score(problem_ab_lowp, [1.0, 1.0, 0.0])
        res_ab_lowp = optimize(problem_ab_lowp)
        println("problem=$problem_ab_lowp")
        @test_skip res_ab_lowp.weights ≈ [0.0, 0.0, 1.0] atol=1E-4

        # higher prior probability to select sets, lower probability to miss active element, so select a and b
        problem_ab = CoverProblem(sm_ab, CoverParams(a=0.1, b=0.1))
        @test score(problem_ab, [1.0, 1.0, 0.0]) < score(problem_ab, [0.0, 0.0, 1.0])
        res_ab = optimize(problem_ab)
        @test res_ab.weights ≈ [1.0, 1.0, 0.0] atol=1E-4
    end
end
