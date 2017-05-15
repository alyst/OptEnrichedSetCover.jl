@testset "logpvalue(A & B)" begin
    using OptEnrichedSetCover.logpvalue

    @test_throws ArgumentError logpvalue(-1, 0, 1, 1)
    @test_throws ArgumentError logpvalue(1, -1, 1, 1, :right)
    @test_throws ArgumentError logpvalue(1, 1, -1, 1)

    @test_throws ArgumentError logpvalue(2, 0, 1, 1)
    @test_throws ArgumentError logpvalue(1, 3, 2, 1, :right)
    @test_throws ArgumentError logpvalue(1, 1, 0, 1)

    @test logpvalue(1, 0, 1, 1, :right) == -Inf
    @test logpvalue(1, 0, 1, 1) == -Inf
    @test logpvalue(1, 0, 1, 1, :left) == 0.0
    @test logpvalue(1, 1, 1, 1, :left) == 0.0
    @test logpvalue(1, 0, 1, 1, :both) == 0.0
    @test logpvalue(1, 1, 1, 1, :both) == 0.0

    @test logpvalue(1, 2, 3, 1, :right) ≈ log(2/3)
    @test logpvalue(1, 2, 3, 1, :left) == 0.0
    @test logpvalue(2, 2, 4, 1, :left) ≈ log(5/6)
    @test logpvalue(2, 2, 4, 1, :right) ≈ log(5/6)

    # no overlap
    @test logpvalue(10, 10, 100, -1, :left) == -Inf
    @test logpvalue(10, 10, 100, 0, :right) == 0.0

    # full overlap
    @test logpvalue(10, 10, 100, 10, :left) == 0.0
    @test logpvalue(10, 10, 100, 10, :right) ≈ -log(binomial(100, 10))
    @test logpvalue(10, 10, 100, 11, :right) == -Inf

    @test logpvalue(10, 10, 20, 5, :left) ≈ log(0.6718591)
    @test logpvalue(10, 10, 20, 5, :right) ≈ log(0.6718591)

    @test logpvalue(15, 15, 30, 10, :right) ≈ log(0.071555489)
    @test logpvalue(15, 15, 30, 10, :left) ≈ log(0.9865811354)

    @test logpvalue(15, 15, 30, 5, :right) ≈ log(0.9865811354)
    @test logpvalue(15, 15, 30, 5, :left) ≈ log(0.071555489)
end
