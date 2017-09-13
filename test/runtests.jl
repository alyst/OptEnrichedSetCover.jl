using OptEnrichedSetCover
using Base.Test, Compat, DataFrames

# write your own tests here

my_tests = [
    "test_set_score.jl",
    "test_sparse_mask_matrix.jl",
    "test_mosaic.jl",
    "test_cover_problem.jl",
    "test_cover_enumerator.jl",
    "test_parallel.jl"
]

for t in my_tests
    info("Testing $t...")
    include(t)
end
