include("parameter.jl")
include("markov.jl")

using .Markov

@assert length(ARGS) >= 1 "Parameters PID, seed are expected!"
PID = parse(Int, ARGS[1])
RNG = length(ARGS) > 1 ? MersenneTwister(parse(Int, ARGS[2])) : MersenneTwister()
# rand(RNG) to generate a random number

printStatus()

println(KGrid.grid)

# ReWeight[2] = 2.0

println(ReWeight)

# println(Para.Counter)

