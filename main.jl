include("parameter.jl")
include("grid.jl")

using Random

@assert length(ARGS) >= 1 "Parameters PID, seed are expected!"
PID = parse(Int, ARGS[1])
RNG = length(ARGS) > 1 ? MersenneTwister(parse(Int, ARGS[2])) : MersenneTwister()

mutable struct State
    step::Int
    rng::MersenneTwister
    order::Int
    scaleidx::Int
    extTidx::Int
    extKidx::Int
    extAngidx::Int
    absWeight::Float
    T::Vector{Float}
    K::Vector{Mom}
    function State(rng)
        varT = rand(rng, LastTidx) .* Beta
        varK = rand(rng, Mom, LastKidx) .* Kf
        # rand!(rng, varK)
        curr = new(0, rng, 0, 1, 1, 1, 1, 0.0, varT, varK)
        curr.T[LastTidx] = Grid.tau.grid[curr.extTidx]
        curr.T[1] = 0.0
        curr.T[2] = Beta / 2.0
        curr.T[3] = Beta / 3.0
        curr.T[4] = Beta / 4.0
        curr.T[5] = Beta / 5.0
        curr.T[6] = Beta / 6.0

        curr.K[5] = [0.0, Kf, 0.0]
        curr.K[6] = [0.0, Kf, 0.0]
        curr.K[7] = [0.0, Kf, 0.0]
        curr.K[8] = [0.0, Kf, 0.0]
        curr.K[9] = [0.0, Kf, 0.0]
        if DiagType == GAMMA
            kL = [Kf, 0.0, 0.0]
            curr.K[OUTL] = curr.K[INL] = kL[1:DIM]
            θ = acos(Grid.angle.grid[curr.extAngidx])
            kR = [Kf * cos(θ), Kf * sin(θ), 0.0]
            curr.K[OUTR] = curr.K[INR] = kR[1:DIM]
        else
            k = [Grid.K.grid[curr.extK], 0.0, 0.0]
            curr.K[1] = k[1:DIM]
        end
        return curr
    end
end

const Curr = State(RNG)
include("markov.jl")
Markov.init()
Markov.test()

import BenchmarkTools: @btime
for _order in 1:Order
    println("Benchmark Order $_order")
    @btime Markov.benchmark(o) samples = 1 evals = 1000 setup = (o = $_order)
    println(sum(Markov.ver4[_order].weight))
end



# Markov.init(Counter, RNG)

Markov.printStatus()
