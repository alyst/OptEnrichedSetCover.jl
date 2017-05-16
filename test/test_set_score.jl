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
    @test logpvalue(1, 0, 1, 1, :both) == -Inf
    @test logpvalue(1, 1, 1, 1, :both) == 0.0
    @test logpvalue(10, 5, 12, 1, :left) == -Inf
    @test logpvalue(10, 5, 12, 1, :right) == 0.0
    @test logpvalue(10, 5, 12, 1, :both) == -Inf
    @test_throws ArgumentError logpvalue(1, 0, 1, 1, :center)

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

    # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=17271
    @test logpvalue(22837, 4, 22837, 4, :left) == 0.0
    @test logpvalue(22837, 4, 22837, 4, :right) == 0.0
    @test logpvalue(22837, 22837-4, 22837, 22837-4, :left) == 0.0
    @test logpvalue(22837, 22837-4, 22837, 22837-4, :right) == 0.0

    @test logpvalue(10, 10, 20, 5, :left) ≈ log(0.6718591)
    @test logpvalue(10, 10, 20, 5, :right) ≈ log(0.6718591)
    @test logpvalue(10, 10, 20, 5, :both) ≈ 0.0

    @test logpvalue(15, 15, 30, 10, :right) ≈ log(0.071555489)
    @test logpvalue(15, 15, 30, 10, :left) ≈ log(0.9865811354)
    @test logpvalue(15, 15, 30, 10, :both) ≈ log(2*0.071555489)

    @test logpvalue(15, 15, 30, 5, :right) ≈ log(0.9865811354)
    @test logpvalue(15, 15, 30, 5, :left) ≈ log(0.071555489)
end
