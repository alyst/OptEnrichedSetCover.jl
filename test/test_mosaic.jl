facts("SetMosaic") do
    context("empty") do
        sm = OESC.SetMosaic(Set{Symbol}[]);

        @fact OESC.nelements(sm) --> 0
        @fact OESC.ntiles(sm) --> 0
        @fact OESC.nsets(sm) --> 0
    end

    context("empty but with :a :b elements") do
        sm = OESC.SetMosaic(Set{Symbol}[], Set([:a, :b]))
        @fact OESC.nelements(sm) --> 2
        @fact OESC.ntiles(sm) --> 0
        @fact OESC.nsets(sm) --> 0
    end

    context("[]") do
        sm = OESC.SetMosaic([Set(Symbol[])]);

        @fact OESC.nelements(sm) --> 0
        @fact OESC.ntiles(sm) --> 0
        @fact OESC.nsets(sm) --> 1
    end

    context("[:a]") do
        sm = OESC.SetMosaic([Set(Symbol[:a])]);

        @fact OESC.nelements(sm) --> 1
        @fact OESC.ntiles(sm) --> 1
        @fact OESC.tile(sm, 1) --> [1]
        @fact OESC.nsets(sm) --> 1
    end

    context("[:a] []") do
        sm = OESC.SetMosaic([Set(Symbol[:a]), Set(Symbol[])]);

        @fact OESC.nelements(sm) --> 1
        @fact OESC.ntiles(sm) --> 1
        @fact OESC.tile(sm, 1) --> [1]
        @fact OESC.nsets(sm) --> 2
    end

    context("[:a] [:b]") do
        sm = OESC.SetMosaic([Set(Symbol[:a]), Set(Symbol[:b])]);

        @fact OESC.nelements(sm) --> 2
        @fact OESC.ntiles(sm) --> 2
        @pending OESC.tile(sm, 1) --> [1] # FIXME
        @pending OESC.tile(sm, 2) --> [2] # FIXME
        @fact OESC.nsets(sm) --> 2
    end

    context("[:a] [:a]") do
        sm = OESC.SetMosaic([Set(Symbol[:a]), Set(Symbol[:a])]);

        @fact OESC.nelements(sm) --> 1
        @fact OESC.ntiles(sm) --> 1
        @fact OESC.tile(sm, 1) --> [1]
        @fact OESC.nsets(sm) --> 2
    end

    context("[:a] [:b] [:a :b]") do
        sm = OESC.SetMosaic([Set(Symbol[:a]), Set(Symbol[:b]), Set(Symbol[:a, :b])]);

        @fact OESC.nelements(sm) --> 2
        @fact OESC.ntiles(sm) --> 2
        @pending OESC.tile(sm, 1) --> [1] # FIXME
        @pending OESC.tile(sm, 2) --> [2] # FIXME
        @fact OESC.nsets(sm) --> 3
    end

    context("[:a :b] [:c :d] [:a :b :c]") do
        sm = OESC.SetMosaic([Set(Symbol[:a, :b]), Set(Symbol[:c, :d]), Set(Symbol[:a, :b, :c])]);

        @fact OESC.nelements(sm) --> 4
        @fact OESC.ntiles(sm) --> 3
        @pending OESC.tile(sm, 1) --> [1, 2] # FIXME
        @pending OESC.tile(sm, 2) --> [3] # FIXME
        @pending OESC.tile(sm, 3) --> [4] # FIXME
        @fact OESC.nsets(sm) --> 3
    end

    context("[:a :b] [:c :d] [:a :b :c :d]") do
        sm = OESC.SetMosaic([Set(Symbol[:a, :b]), Set(Symbol[:c, :d]), Set(Symbol[:a, :b, :c, :d])]);

        @fact OESC.nelements(sm) --> 4
        @fact OESC.ntiles(sm) --> 2
        @pending OESC.tile(sm, 1) --> [1, 2] # FIXME
        @pending OESC.tile(sm, 2) --> [3, 4] # FIXME
        @fact OESC.nsets(sm) --> 3
    end
end

facts("MaskedSetMosaic") do
    context("empty") do
        sm = OESC.SetMosaic(Set{Symbol}[])
        msm = OESC.mask(sm, Set{Symbol}())
        @fact OESC.unmask(msm) --> sm
        @fact OESC.nelements(msm) --> 0
        @fact OESC.ntiles(msm) --> 0
        @fact OESC.nsets(msm) --> 0
        @fact OESC.nmasked_pertile(msm) --> Int[]
        @fact OESC.nunmasked_pertile(msm) --> Int[]
    end

    context("empty but with elements") do
        sm = OESC.SetMosaic(Set{Symbol}[], Set([:a, :b]))
        msm = OESC.mask(sm, Set([:a]))
        @fact OESC.unmask(msm) --> sm
        @fact OESC.nelements(msm) --> 2
        @fact OESC.ntiles(msm) --> 0
        @fact OESC.nsets(msm) --> 0
        @fact OESC.nmasked(msm) --> 1
        @fact OESC.nunmasked(msm) --> 1
        @fact OESC.nmasked_pertile(msm) --> Int[]
        @fact OESC.nunmasked_pertile(msm) --> Int[]
    end

    context("[:a :b] [:c :d] [:a :b :c :d], :a :b") do
        sm = OESC.SetMosaic([Set(Symbol[:a, :b]), Set(Symbol[:c, :d]), Set(Symbol[:a, :b, :c, :d])]);

        msm = OESC.mask(sm, Set(Symbol[:a, :b]))

        @fact OESC.nelements(msm) --> 4
        @fact OESC.ntiles(msm) --> 2
        @fact OESC.nmasked(msm) --> 2
        @fact OESC.nunmasked(msm) --> 2
        @pending OESC.tiles(msm) --> Vector{Int}[Int[1, 2], Int[3, 4]]
        @fact OESC.nsets(msm) --> 2
        @fact length(OESC.nmasked_pertile(msm)) --> 2
        @fact length(OESC.nmasked_perset(msm)) --> 2
        # FIXME set indices not stable
        @pending OESC.nmasked_pertile(msm) --> Int[2, 0]
        @pending OESC.nunmasked_pertile(msm) --> Int[0, 2]
        @pending OESC.nmasked_perset(msm) --> Int[2, 2]
        @pending OESC.nunmasked_perset(msm) --> Int[0, 2]

        msm_copy = copy(msm)
        @fact OESC.nelements(msm_copy) --> OESC.nelements(msm)
        @fact OESC.nmasked(msm_copy) --> OESC.nmasked(msm)
        @fact msm_copy.original --> msm.original
        @fact OESC.nsets(msm_copy) --> OESC.nsets(msm)

        # mask with nonexisting element
        msm2 = OESC.mask(sm, Set(Symbol[:a, :b, :g]))
        @fact OESC.nelements(msm2) --> 4
        @fact OESC.ntiles(msm2) --> 2
        @fact OESC.nmasked(msm2) --> 2
        @fact OESC.nunmasked(msm2) --> 2
    end

    context("A=[:a :b] B=[:c :d] C=[:a :b :c :d], :a :b") do
        sm = OESC.SetMosaic(Dict(:A=>Set(Symbol[:a, :b]), :B=>Set(Symbol[:c, :d]), :C=>Set(Symbol[:a, :b, :c, :d])));
        msm = OESC.mask(sm, Set(Symbol[:a, :b]))

        @fact OESC.nelements(msm) --> 4
        @fact OESC.ntiles(msm) --> 2
        @pending OESC.tiles(msm) --> Vector{Int}[Int[1, 2], Int[3, 4]]
        @fact OESC.nsets(msm) --> 2
        @pending OESC.nmasked_pertile(msm) --> Int[2, 0]
        @pending OESC.nunmasked_pertile(msm) --> Int[0, 2]
    end

    context("filter!()") do
        sm = OESC.SetMosaic([Set(Symbol[:a, :b]), Set(Symbol[:c, :d]), Set(Symbol[:a, :b, :c, :d])]);
        msm = OESC.mask(sm, Set(Symbol[:b, :c]))

        @fact OESC.nelements(msm) --> 4
        @fact OESC.ntiles(msm) --> 2
        @fact OESC.nsets(msm) --> 3
        @fact OESC.nmasked_pertile(msm) --> Int[1, 1]
        @fact OESC.nunmasked_pertile(msm) --> Int[1, 1]

        filter!(msm, Bool[true, true, false])
        @fact OESC.nelements(msm) --> 4
        @fact OESC.ntiles(msm) --> 2
        @fact OESC.nsets(msm) --> 2
        @fact OESC.nmasked_pertile(msm) --> Int[1, 1]
        @fact OESC.nunmasked_pertile(msm) --> Int[1, 1]

        filter!(msm, Bool[false, true])
        @fact OESC.nelements(msm) --> 4
        @fact OESC.ntiles(msm) --> 1
        @fact OESC.nsets(msm) --> 1
        @fact OESC.nmasked_pertile(msm) --> Int[1]
        @fact OESC.nunmasked_pertile(msm) --> Int[1]
    end
end
