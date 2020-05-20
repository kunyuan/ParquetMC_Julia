include("parameter.jl")

PID = parse(Int, ARGS[1])
if length(ARGS) > 1
    RNG = MersenneTwister(parse(Int, ARGS[2])); # rand(RNG) to generate a random number
end
# Random.seed!(seed)
println(rand(RNG))

println(KGrid.grid)
println(Counter)
