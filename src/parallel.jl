immutable CollectCovers{MC,SC}
    mosaics::MC
    sets::SC
    cover_params::CoverParams
    max_set_score::Float64

    CollectCovers(worker::CollectCovers{MC,SC}) =
        new(worker.mosaics, worker.sets,
            worker.cover_params, worker.max_set_score)

    function (::Type{CollectCovers}){MC,SC}(
        mosaics::MC, sets::SC,
        cover_params::CoverParams=CoverParams(sel_prob=0.1),
        max_set_score::Number=log(1E-3)
    )
        new{MC,SC}(mosaics, sets, cover_params, max_set_score)
    end
end

function (worker::CollectCovers)(mosaic_key, set_key)
    mosaic_masked = mask(worker.mosaics[mosaic_key], worker.sets[set_key])
    covers_enum = CoverEnumerator(mosaic_masked, worker.cover_params)
    covers_coll = collect(covers_enum; max_set_score=worker.max_set_score)
    ((mosaic_key, set_key), !isempty(covers_coll) ? covers_coll : nothing)
end

function pcollect{MK,SK}(
    mosaics::Dict{MK}, sets::Dict{SK};
    cover_params::CoverParams=CoverParams(sel_prob=0.1),
    max_set_score::Float64=log(1E-3), pids=workers(), mode=:parallel,
    verbose::Bool=false
)
    info("Parallel OESC ($(length(mosaics)) mosaic(s) × $(length(sets)) set(s))...")
    # FIXME workaround for pmap() being unable to handle callable objects +
    # unable to efficiently handle large callable objects
    collect_covers = CollectCovers(mosaics, sets, cover_params, max_set_score)
    #raw_res = sizehint!(Vector{Tuple{MK,SK,CoverCollection}}(),
    #                    length(tasks))
    mosaicXset_keys = vec([(m, s) for m in keys(mosaics), s in keys(sets)])
    if mode == :parallel
        keyXcover_vec = @parallel (append!) for (m,s) in shuffle!(mosaicXset_keys)
            res = Any[collect_covers(m, s)]
            info("$(myid()): done mosaic=$m × set=$s")
            res
        end
    elseif mode == :sequential
        keyXcover_vec = map(mosaicXset_keys) do mXs
            m, s = mXs
            verbose && info("Calculating mosaic=$m × set=$s...")
            collect_covers(m, s)
        end
    else
        throw(ArgumentError("Unknown mode $mode"))
    end
    verbose && (@show raw_res)

    verbose && info("Done covers enumeration, preparing the results...")
    # filter out empty covers
    keyXcover = sizehint!(Dict{Tuple{MK,SK}, CoverCollection}(), length(keyXcover_vec))
    for (mXs_key, covers_coll) in keyXcover_vec
        if covers_coll !== nothing
            keyXcover[mXs_key] = covers_coll
        end
    end
    info("$(length(keyXcover)) non-trivial cover(s) collected")
    return keyXcover
end

function pcollect(sets1_colls, sets2;
    max_set_score::Float64=log(1E-3),
    cover_params::CoverParams=CoverParams(sel_prob=0.1),
    mode::Symbol=:sequential, verbose::Bool=false
)
    verbose && info("Preparing mosaics...")
    #p = Progress(sum(length, values(entrez_clusters))*length(entrez_colls), 5.0,
    #              "Computing cluster coverings...")   # minimum update interval: 5 second
    sets1_mosaics = Dict(pmap(sets1_colls) do coll_kv
        verbose && info("Preparing $(coll_kv[1]) mosaic...")
        Pair(coll_kv[1], SetMosaic(coll_kv[2]))
    end)
    return pcollect(sets1_mosaics, sets2,
                    max_set_score=max_set_score,
                    cover_params=cover_params, mode=mode,
                    verbose=verbose)
end
