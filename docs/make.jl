using OptEnrichedSetCover, DataFrames
using Documenter

makedocs(
    format = Documenter.HTML(prettyurls=false),
    sitename = "OptEnrichedSetCover.jl",
    authors = "Alexey Stukalov",
    modules = [OptEnrichedSetCover],
    clean = true,
    pages = [
        "Introduction" => "index.md",
        "Method description" => "method.md",
        "Sets collection" => "mosaic.md",
        "Cover problem" => "cover_problem.md",
        "Internals" => [
            "Array pool" => "array_pool.md",
            "Sparse mask matrix" => "sparse_mask_matrix.md",
        ]
    ],
)

deploydocs(
    repo = "github.com/alyst/OptEnrichedSetCover.jl.git",
)
