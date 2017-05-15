@testset "logpvalue(A & B)" begin
    using OptEnrichedSetCover.logpvalue

    @test_throws ArgumentError logpvalue(-1, 0, 1, 1)
    @test_throws ArgumentError logpvalue(1, -1, 1, 1, :left)
    @test_throws ArgumentError logpvalue(1, 1, -1, 1)

    @test_throws ArgumentError logpvalue(2, 0, 1, 1)
    @test_throws ArgumentError logpvalue(1, 3, 2, 1, :left)
    @test_throws ArgumentError logpvalue(1, 1, 0, 1)

    @test_broken logpvalue(1, 0, 1, 1, :left) == -Inf
    @test_broken logpvalue(1, 0, 1, 1) == -Inf
    @test_broken logpvalue(1, 0, 1, 1, :right) == 0.0
    @test_broken logpvalue(1, 1, 1, 1, :right) == 0.0
    @test logpvalue(1, 0, 1, 1, :both) == 0.0
    @test logpvalue(1, 1, 1, 1, :both) == 0.0

    @test_broken logpvalue(1, 2, 3, 1, :left) ≈ log(2/3)
    @test_broken logpvalue(1, 2, 3, 1, :right) == 0.0
    @test logpvalue(2, 2, 4, 1, :left) ≈ log(5/6)
    @test logpvalue(2, 2, 4, 1, :right) ≈ log(5/6)

    @test logpvalue(10, 10, 20, 5, :left) ≈ log(0.6718591)
    @test logpvalue(10, 10, 20, 5, :right) ≈ log(0.6718591)

    @test logpvalue(15, 15, 30, 10, :left) ≈ log(0.071555489)
    @test logpvalue(15, 15, 30, 10, :right) ≈ log(0.9865811354)

    @test logpvalue(15, 15, 30, 5, :left) ≈ log(0.9865811354)
    @test logpvalue(15, 15, 30, 5, :right) ≈ log(0.071555489)
end
