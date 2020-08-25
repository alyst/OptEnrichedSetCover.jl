# [Array pool](@id array_pool)

The `ArrayPool` type is meant to reduce the stress on Julia *Garbage Collector*
by maintaining the pool of preallocated arrays that the program can temporarily
borrow and then return back to the pool when not needed anymore.

```@autodocs
Modules=[OptEnrichedSetCover]
Pages = ["array_pool.jl"]
```
