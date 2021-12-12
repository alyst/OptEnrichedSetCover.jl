using OptEnrichedSetCover
using Test, LinearAlgebra, DataFrames

const OESC = OptEnrichedSetCover

# run the tests from jlfile
macro testfile(jlfile)
    quote
        @testset "\"$($jlfile)\" tests" begin
            include($jlfile)
        end
    end
end

@testfile "test_set_score.jl"
@testfile "test_sparse_mask_matrix.jl"
@testfile "test_mosaic.jl"
@testfile "test_masked_mosaic.jl"
@testfile "test_weighted_mosaic.jl"
@testfile "test_multiobj_cover_problem.jl"
@testfile "test_cover_enumerator.jl"
@testfile "test_parallel.jl"
