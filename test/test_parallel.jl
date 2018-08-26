global pcollect_opt_methods = [MultiobjectiveOptimizerParams()]
if OptEnrichedSetCover.__quadratic_problem_supported__
    push!(pcollect_opt_methods, QuadraticOptimizerParams())
else
    @warn "pcollect(): Quadratic problem not supported and not tested"
end

@testset "pcollect(opt_params=$opt_params)" for opt_params in pcollect_opt_methods
    mosaics = Dict("A" => SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                                    Set([:a, :b, :c, :d, :e])))
    masks = Dict(1 => Set([:a, :b]),
                 2 => Set([:d, :e]),
                 3 => Set([:a, :c, :d]))
    res = OptEnrichedSetCover.pcollect(mosaics, masks, mode=:sequential,
                                       opt_params=opt_params,
                                       cover_params=CoverParams(sel_prob=0.75),
                                       enum_params=CoverEnumerationParams(max_set_score=10.0))
    @test length(res) == 2
    @test isa(res, Dict{Tuple{String, Int}, CoverCollection})
end
