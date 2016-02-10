facts("logpvalue(A & B)") do
    @fact_throws OESC.logpvalue(-1, 0, 1, 1, tail=:left)
    @fact_throws OESC.logpvalue(1, -1, 1, 1, tail=:left)
    @fact_throws OESC.logpvalue(1, 1, -1, 1, tail=:left)

    @fact_throws OESC.logpvalue(2, 0, 1, 1, tail=:left)
    @fact_throws OESC.logpvalue(1, 3, 2, 1, tail=:left)
    @fact_throws OESC.logpvalue(1, 1, 0, 1, tail=:left)

    @fact OESC.logpvalue(1, 0, 1, 1, tail=:left) --> -Inf
    @fact OESC.logpvalue(1, 0, 1, 1, tail=:right) --> 0.0
    @fact OESC.logpvalue(1, 1, 1, 1, tail=:right) --> 0.0

    @fact OESC.logpvalue(1, 2, 3, 1, tail=:left) --> roughly(log(2/3))
    @fact OESC.logpvalue(1, 2, 3, 1, tail=:right) --> roughly(0.0)
    @fact OESC.logpvalue(2, 2, 4, 1, tail=:left) --> roughly(log(5/6))
    @fact OESC.logpvalue(2, 2, 4, 1, tail=:right) --> roughly(log(5/6))

    @fact OESC.logpvalue(10, 10, 20, 5, tail=:left) --> roughly(log(0.6718591))
    @fact OESC.logpvalue(10, 10, 20, 5, tail=:right) --> roughly(log(0.6718591))

    @fact OESC.logpvalue(15, 15, 30, 10, tail=:left) --> roughly(log(0.071555489))
    @fact OESC.logpvalue(15, 15, 30, 10, tail=:right) --> roughly(log(0.9865811354))

    @fact OESC.logpvalue(15, 15, 30, 5, tail=:left) --> roughly(log(0.9865811354))
    @fact OESC.logpvalue(15, 15, 30, 5, tail=:right) --> roughly(log(0.071555489))
end
