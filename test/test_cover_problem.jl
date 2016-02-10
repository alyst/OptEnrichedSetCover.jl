facts("CoverProblem") do
    context("CoverParams") do
    end

    context("empty") do
        sm = OESC.SetMosaic(Set{Symbol}[])
        msm = OESC.mask(sm, Set{Symbol}())
        problem = OESC.CoverProblem(msm, OESC.CoverParams(0.1, 0.2, 0.5))
        @fact OESC.nsets(problem) --> 0
        @fact OESC.score(problem, Float64[]) --> 0.0

        res = OESC.optimize(problem)
        @fact res.weights --> Float64[]
        @fact res.score --> 0.0
    end

    context("[:a]") do # FIXME take element detection probability into account
        # empty mask problem
        empty_problem = OESC.CoverProblem(OESC.mask(OESC.SetMosaic([Set([:a])]), Set{Symbol}()), OESC.CoverParams(0.5, 0.5))
        @fact OESC.nsets(empty_problem) --> 0
        empty_res = OESC.optimize(empty_problem)
        @fact empty_res.weights --> Float64[]
        @fact empty_res.score --> 0.0

        # disabled because :a is all elements
        dis_problem = OESC.CoverProblem(OESC.mask(OESC.SetMosaic([Set([:a])], Set([:a])), Set{Symbol}([:a])),
                                        OESC.CoverParams(0.5, 0.5))
        @fact OESC.nsets(dis_problem) --> 1
        dis_res = OESC.optimize(dis_problem)
        @fact dis_res.weights --> roughly(Float64[0.0], atol=1E-8)

        # enabled because there is :a and :b
        en_problem = OESC.CoverProblem(OESC.mask(OESC.SetMosaic([Set([:a])], Set([:a, :b])), Set{Symbol}([:a])),
                                       OESC.CoverParams(0.5, 0.5))
        @fact OESC.nsets(en_problem) --> 1
        en_res = OESC.optimize(en_problem)
        @fact en_res.weights --> roughly(Float64[1.0], atol=1E-8)
    end

    context("[:a :b] [:c :d] [:a :b :c :d], mask=[:a :b]") do # FIXME take weights into account
        sm = OESC.SetMosaic([Set([:a, :b]), Set([:c, :d]), Set([:a, :b, :c, :d])]);
        sm_ab = OESC.mask(sm, Set(Symbol[:a, :b]))

        problem_ab_lowp = OESC.CoverProblem(sm_ab, OESC.CoverParams(0.1, 0.1))
        @pending OESC.nsets(problem_ab_lowp) --> 2

        res_ab_lowp = OESC.optimize(problem_ab_lowp)
        @pending res_ab_lowp.weights --> roughly(Float64[0.0, 0.0], atol=1E-4)

        problem_ab = OESC.CoverProblem(sm_ab, OESC.CoverParams(0.1, 0.1))
        res_ab = OESC.optimize(problem_ab)
        @fact res_ab.weights --> roughly(Float64[1.0, 0.0], atol=1E-4)
    end

    context("[:a :b] [:b :c] [:a :b :c], mask=[:b]") do # FIXME take weights into account
        sm = OESC.SetMosaic([Set([:a, :b]), Set([:b, :c]), Set([:a, :b, :c])])
        sm_b = OESC.mask(sm, Set([:b]))

        problem_b_lowp = OESC.CoverProblem(sm_b, OESC.CoverParams(0.1, 0.1))
        @fact OESC.nsets(problem_b_lowp) --> 3

        res_b_lowp = OESC.optimize(problem_b_lowp)
        @pending res_b_lowp.weights --> roughly(Float64[0.0, 0.0, 0.0], atol=1E-4)

        problem_b = OESC.CoverProblem(sm_b, OESC.CoverParams(0.4, 0.1))
        res_b = OESC.optimize(problem_b)
        # [:a :b] is as good as [:b :c]
        @fact res_b.weights[1] --> roughly(res_b.weights[2])
        @fact res_b.weights[3] --> roughly(0.0)
    end

    context("[:a] [:b] [::c] [:a :b :c], mask=[:a :b]") do # FIXME take weights into account
        sm = OESC.SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                            Set([:a, :b, :c, :d, :e]))
        sm_ab = OESC.mask(sm, Set([:a, :b]))

        # low prior probability to select sets, high probability to miss active element, so select abc
        problem_ab_lowp = OESC.CoverProblem(sm_ab, OESC.CoverParams(0.4, 0.1))
        @fact OESC.nsets(problem_ab_lowp) --> 3 # c is out
        @pending OESC.score(problem_ab_lowp, [0.0, 0.0, 1.0]) --> less_than(OESC.score(problem_ab_lowp, [1.0, 1.0, 0.0]))
        res_ab_lowp = OESC.optimize(problem_ab_lowp)
        println("problem=$problem_ab_lowp")
        @pending res_ab_lowp.weights --> roughly(Float64[0.0, 0.0, 1.0], atol=1E-4)

        # higher prior probability to select sets, lower probability to miss active element, so select a and b
        problem_ab = OESC.CoverProblem(sm_ab, OESC.CoverParams(0.1, 0.1))
        @fact OESC.score(problem_ab, [1.0, 1.0, 0.0]) --> less_than(OESC.score(problem_ab, [0.0, 0.0, 1.0]))
        res_ab = OESC.optimize(problem_ab)
        @fact res_ab.weights --> roughly(Float64[1.0, 1.0, 0.0], atol=1E-4)
    end
end
