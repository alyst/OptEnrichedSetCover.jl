include("/home/astukalov/projects/oesc_jl/test/runtests.jl")

sets = [Set([:a]), Set([:b]), Set([:a :b]), Set([:a, :b, :c]), Set([:d, :e, :f]), Set([:a, :b, :c, :d, :e, :f])]
active = Set([:a, :b, :d])

include("/home/astukalov/projects/oesc_jl/oesc.jl")
params = OESC.CoverParams(0.1, 0.01, 0.01);
sm = OESC.SetMosaic(sets);
msm = OESC.mask(sm, Set([:a,:b,:e,:c]));
problem = OESC.CoverProblem(msm, params);
opt_model = OESC.opt_model(problem)

using JuMP
solve(opt_model)
getValue(getVar(opt_model, :w))

res = OESC.optimize(problem)

using Gadfly
plot(x=linspace(0.0,1.0,100), y=linspace(0.0,1.0,100),
     z=(x,y) -> OESC.score(problem, 1.0 - [x,y]),
     Geom.contour(levels=50))

OESC.score(problem, res.set_weights)

sm = OESC.SetMosaic(Set{Symbol}[Set(Symbol[:a])])
msm = OESC.mask(sm, Set{Symbol}([:a]))

en_problem = OESC.CoverProblem(msm, OESC.CoverParams(0.5, 0.5, 0.7));
plot(x -> OESC.score(en_problem, 1.0 - [x]), 0.0, 1.0)

dis_problem = OESC.CoverProblem(msm, OESC.CoverParams(0.5, 0.5, 0.4));
plot(x -> OESC.score(dis_problem, 1.0 - [x]), 0.0, 1.0)

en_res = OESC.optimize(en_problem)
@fact en_res.set_weights --> roughly(Float64[1.0], atol=1E-8)

sm = OESC.SetMosaic([Set(Symbol[:a, :b]), Set(Symbol[:b, :c]), Set(Symbol[:a, :b, :c])])
sm_ab = OESC.mask(sm, Set(Symbol[:a]))

problem_ab_lowp = OESC.CoverProblem(sm_ab, OESC.CoverParams(0.1, 0.1, 0.01))
@fact OESC.nsets(problem_ab_lowp) --> 3
@fact OESC.ntiles(problem_ab_lowp) --> 3

res_ab_lowp = OESC.optimize(problem_ab_lowp)
@fact res_ab_lowp.set_weights --> roughly(Float64[0.0, 0.0, 0.0], atol=1E-4)

problem_ab = OESC.CoverProblem(sm_ab, OESC.CoverParams(0.4, 0.1, 0.3))
res_ab = OESC.optimize(problem_ab)
@fact res_ab.set_weights --> roughly(Float64[1.0, 0.0, 0.0], atol=1E-4)

sm = OESC.SetMosaic([Set(Symbol[:a]), Set(Symbol[:b]),
                      Set(Symbol[:c]), Set(Symbol[:a, :b, :c])])
sm_ab = OESC.mask(sm, Set(Symbol[:a, :b]))

# low prior probability to select sets, high probability to miss active element, so select abc
problem_ab = OESC.CoverProblem(sm_ab, OESC.CoverParams(0.5, 0.1, 0.1))
OESC.score(problem_ab, [1.0, 1.0, 0.0])

using Gadfly
plot(x=linspace(0.0,1.0,100), y=linspace(0.0,1.0,100),
     z=(x,y) -> OESC.score(problem_ab, Float64[x,y,0.001]),
     Geom.contour(levels=50))

include("/home/astukalov/projects/oesc/oesc.jl")
sm = OESC.SetMosaic([Set(Symbol[:a]), Set(Symbol[:b]),
                      Set(Symbol[:c]), Set(Symbol[:a, :b, :c])])
sm_ab = OESC.mask(sm, Set(Symbol[:a, :b]))

cover_etor = OESC.CoverEnumerator(sm_ab, OESC.CoverParams(0.5, 0.1, 0.1))
cover_coll = collect(cover_etor)
DataFrame(cover_coll)

cover_etor2 = OESC.CoverEnumerator(sm_ab, OESC.CoverParams(0.1, 0.1, 0.2))
cover_coll2 = collect(cover_etor2)

cover_results[(:MF, 1, :-, 3796)]

@where(cover_results_df, :prob .> 0.8)

ProfileView.view()

clu_go_bp = OESC.mask(go_bp, exp_comp_clusters["2+"][7])
OESC.nsets(clu_go_bp)
OESC.ntiles(clu_go_bp)

cover_etor = OESC.CoverEnumerator(clu_go_bp, OESC.CoverParams(0.5, 0.1, 0.1))
cover_coll = collect(cover_etor)
