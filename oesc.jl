# Optimal Enriched-Set Cover (OESC)

module OESC

using DataFrames, Distributions, MathProgBase, JuMP, Ipopt

for script_file in [ "set_score.jl",
                     "sparse_mask_matrix.jl",
                     "mosaic.jl",
                     "cover_problem.jl",
                     "cover_enumerator.jl" ]
    include(script_file)
end

end
