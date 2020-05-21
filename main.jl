include("parameter.jl")
include("markov.jl")

using Random

# using .Markov

@assert length(ARGS) >= 1 "Parameters PID, seed are expected!"
PID = parse(Int, ARGS[1])
RNG = length(ARGS) > 1 ? MersenneTwister(parse(Int, ARGS[2])) : MersenneTwister()

println(rand(RNG))

counter = Counter()
counter.step()
println(counter.state())
Markov.init(counter, RNG)
println(counter.state())


# Markov.init(Counter, RNG)

Markov.printStatus()
