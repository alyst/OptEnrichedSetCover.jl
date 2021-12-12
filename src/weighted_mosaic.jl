"""
[`SetMosaic`](@ref) with the weights assigned to its sets.

## Type parameters
* `T`: the type of elements
* `S`: the type of set ids
* `E`: the type of experiment ids
* `W`: type of the weight
"""
abstract type AbstractWeightedSetMosaic{T,S,E,W} end

"""
    originalmosaic(mosaic::AbstractWeightedSetMosaic) -> SetMosaic

Get the original [`SetMosaic`](@ref).
"""
originalmosaic(mosaic::AbstractWeightedSetMosaic) = mosaic.original

Base.eltype(::Type{<:AbstractWeightedSetMosaic{T}}) where T = T
Base.eltype(mosaic::AbstractWeightedSetMosaic) = eltype(typeof(mosaic))

setidtype(::Type{<:AbstractWeightedSetMosaic{<:Any, S}}) where S = S
setidtype(mosaic::AbstractWeightedSetMosaic) = setidtype(typeof(mosaic))

experimentidtype(::Type{<:AbstractWeightedSetMosaic{<:Any, <:Any, E}}) where E = E
experimentidtype(mosaic::AbstractWeightedSetMosaic) = experimentidtype(typeof(mosaic))

weighttype(::Type{<:AbstractWeightedSetMosaic{<:Any, <:Any, <:Any, W}}) where W = W
weighttype(mosaic::AbstractWeightedSetMosaic) = weighttype(typeof(mosaic))
