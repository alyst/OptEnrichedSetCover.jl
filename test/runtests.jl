using OptEnrichedSetCover
using Test, LinearAlgebra, DataFrames

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
if OptEnrichedSetCover.__quadratic_problem_supported__
    @testfile "test_cover_problem.jl"
else
    @warn "test_cover_problem.jl tests skipped"
end
@testfile "test_multiobj_cover_problem.jl"
@testfile "test_cover_enumerator.jl"
@testfile "test_parallel.jl"
