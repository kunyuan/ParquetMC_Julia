include("parameter.jl")
include("markov.jl")

using Random

# using .Markov

@assert length(ARGS) >= 1 "Parameters PID, seed are expected!"
PID = parse(Int, ARGS[1])
RNG = length(ARGS) > 1 ? MersenneTwister(parse(Int, ARGS[2])) : MersenneTwister()
# rand(RNG) to generate a random number


# mutable struct 
#     pid::Int
#     rng::MersenneTwister
#     counter::Int128
#     function Var()
#         pid = parse(Int, ARGS[1])
#         rng = length(ARGS) > 1 ? MersenneTwister(parse(Int, ARGS[2])) : MersenneTwister()
#         counter = 0
#         var = new(pid, rng, counter)
#         return var
#     end
# end

# var = Var()
# println(var.counter)
# var.counter = 2
# println(var.counter)
# Counter = Int128(0)

println(rand(RNG))

counter = Counter()
counter.step()
println(counter.state())
Markov.init(counter, RNG)
println(counter.state())


# Markov.init(Counter, RNG)

Markov.printStatus()

# printStatus()

# println(KGrid.grid)

# ReWeight[2] = 2.0

# println(ReWeight)

# println(Para.Counter)
