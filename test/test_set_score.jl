@testset "logpvalue(A & B)" begin
    using OptEnrichedSetCover.logpvalue

    @test_throws ArgumentError logpvalue(1, -1, 0, 1)
    @test_throws ArgumentError logpvalue(1, 1, -1, 1, :right)
    @test_throws ArgumentError logpvalue(1, 1, 1, -1)

    @test_throws ArgumentError logpvalue(1, 2, 0, 1)
    @test_throws ArgumentError logpvalue(1, 1, 3, 2, :right)
    @test_throws ArgumentError logpvalue(1, 1, 1, 0)

    @test logpvalue(1, 1, 0, 1, :right) == -Inf
    @test logpvalue(1, 1, 0, 1) == -Inf
    @test logpvalue(1, 1, 0, 1, :left) == 0.0
    @test logpvalue(1, 1, 1, 1, :left) == 0.0
    @test logpvalue(1, 1, 0, 1, :both) == -Inf
    @test logpvalue(1, 1, 1, 1, :both) == 0.0
    @test logpvalue(1, 10, 5, 12, :left) == -Inf
    @test logpvalue(1, 10, 5, 12, :right) == 0.0
    @test logpvalue(1, 10, 5, 12, :both) == -Inf
    @test_throws ArgumentError logpvalue(1, 1, 0, 1, :center)

    @test logpvalue(1, 1, 2, 3, :right) ≈ log(2/3)
    @test logpvalue(1, 1, 2, 3, :left) == 0.0
    @test logpvalue(1, 2, 2, 4, :left) ≈ log(5/6)
    @test logpvalue(1, 2, 2, 4, :right) ≈ log(5/6)

    # no overlap
    @test logpvalue(-1, 10, 10, 100, :left) == -Inf
    @test logpvalue( 0, 10, 10, 100, :right) == 0.0

    # full overlap
    @test logpvalue(10, 10, 10, 100, :left) == 0.0
    @test logpvalue(10, 10, 10, 100, :right) ≈ -log(binomial(100, 10))
    @test logpvalue(11, 10, 10, 100, :right) == -Inf

    # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=17271
    @test logpvalue(4, 22837, 4, 22837, :left) == 0.0
    @test logpvalue(4, 22837, 4, 22837, :right) == 0.0
    @test logpvalue(22837-4, 22837, 22837-4, 22837, :left) == 0.0
    @test logpvalue(22837-4, 22837, 22837-4, 22837, :right) == 0.0

    @test logpvalue(5, 10, 10, 20, :left) ≈ log(0.6718591)
    @test logpvalue(5, 10, 10, 20, :right) ≈ log(0.6718591)
    @test logpvalue(5, 10, 10, 20, :both) ≈ 0.0

    @test logpvalue(10, 15, 15, 30, :right) ≈ log(0.071555489)
    @test logpvalue(10, 15, 15, 30, :left) ≈ log(0.9865811354)
    @test logpvalue(10, 15, 15, 30, :both) ≈ log(2*0.071555489)

    @test logpvalue(5, 15, 15, 30, :right) ≈ log(0.9865811354)
    @test logpvalue(5, 15, 15, 30, :left) ≈ log(0.071555489)
end
