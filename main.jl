include("parameter.jl")
include("markov.jl")

using Random

# using .Markov

@assert length(ARGS) >= 1 "Parameters PID, seed are expected!"
PID = parse(Int, ARGS[1])
RNG = length(ARGS) > 1 ? MersenneTwister(parse(Int, ARGS[2])) : MersenneTwister()
println(Main.step)
counter = Counter()
counter.step()

Markov.init()

Markov.test()


# Markov.init(Counter, RNG)

Markov.printStatus()
