@testset "pcollect()" begin
    mosaics = Dict("A" => SetMosaic([Set([:a]), Set([:b]), Set([:c]), Set([:a, :b, :c])],
                                    Set([:a, :b, :c, :d, :e])))
    masks = Dict(1 => Set([:a, :b]),
                 2 => Set([:d, :e]),
                 3 => Set([:a, :c, :d]))
    res = OptEnrichedSetCover.pcollect(mosaics, masks, mode=:sequential, max_set_score=0.0)
    @test length(res) == 2
    @test isa(res, Dict{Tuple{String, Int}, CoverCollection})
end
