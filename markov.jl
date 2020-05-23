module Markov
include("parameter.jl")
include("grid.jl")
include("diag/vertex4.jl")
include("diag/polar.jl")
include("diag/vertex4_test.jl")
using Random, StaticArrays
import .Vertex4, .Ver4Test, .Polar
const UpdateNum = 9

mutable struct State
    step::Int
    order::Int
    absWeight::Float
    scaleidx::Int
    extTidx::Int
    extKidx::Int
    extAngidx::Int
    State() = new(0, 0, 0.0, 1, 1, 1, 1)
end

const curr = State()
const varT = Vector{Float}(undef, LastTidx)
const varK = Vector{Mom}(undef, LastKidx)
const vertex4 = Vector{Vertex4.Ver4}(undef, 0)

# function init(_counter, _rng)
function init()

    # global rng = _rng
    # global counter = _counter

    rng = MersenneTwister()
    counter = Counter()

    #######  initialize MC variables  ################################

    varT .= rand(rng, LastTidx) * Beta
    varT[LastTidx] = Grid.tau.grid[curr.extTidx]
    varT[1] = 0.0
    varT[2] = Beta / 2.0
    varT[3] = Beta / 3.0
    varT[4] = Beta / 4.0
    varT[5] = Beta / 5.0
    varT[6] = Beta / 6.0

    rand!(rng, varK)
    varK .*= Kf
    k0 = varK[1]

    varK[5] = [0.0, Kf, 0.0]
    varK[6] = [0.0, Kf, 0.0]
    varK[7] = [0.0, Kf, 0.0]
    varK[8] = [0.0, Kf, 0.0]
    varK[9] = [0.0, Kf, 0.0]

    if DiagType == GAMMA
        kL = [Kf, 0.0, 0.0]
        varK[OUTL] = varK[INL] = kL[1:DIM]
        θ = acos(Grid.angle.grid[curr.extAngidx])
        kR = [Kf * cos(θ), Kf * sin(θ), 0.0]
        varK[OUTR] = varK[INR] = kR[1:DIM]
    else
        k = [Grid.K.grid[curr.extKidx], 0.0, 0.0]
        varK[1] = k[1:DIM]
    end

    ###### initialized diagram trees #######################################
    # Vertex4.init(varT, varK)
    # Polar.init(varT, varK)
    # Ver4Test.init(varT, varK)

    if DiagType == GAMMA
        chan = [I, T, U, S, TC, UC]
        # chan = [T, ]
        for order = 1:Order
            push!(vertex4, Vertex4.Ver4(0, order, chan, 1, RIGHT, false, true))
        end
        # Vertex4.visualize(vertex4[2])

        GC.gc() # collect garabage, reduce memory usage

        for order in 1:Order
            println("Order $order")
            for i in 1:10
                Vertex4.eval(vertex4[order], varK[INL], varK[OUTL], varK[INR], varK[OUTR], 5, varT, varK, true)
                # println(sum(vertex4[order].weight))
                # println((vertex4[order].weight))
            end
        end

        for order in 1:Order
            @time begin
            # @allocated begin
                println("Order $order")
                for i in 1:1000
                    Vertex4.eval(vertex4[order], varK[INL], varK[OUTL], varK[INR], varK[OUTR], 5, varT, varK, true)
                end
            end
        end
    elseif DiagType == POLAR
        global polar = Vector{Polar.Polarization}(undef, 0)
        for order = 0:Order
            push!(polar, Polar.Polarization(order))
        end
    else
        throw("Not implemented!")
    end
end

function test()
    if DiagType == GAMMA
        Ver4Test.testOneLoopVer4(varT, varK)
    end
end

function printStatus()
    # Var.counter += 1
    ReWeight[2] = 1.6
    println(curr.extAngidx)
    println(Beta)
end

function eval()
    return
end

function changeOrder()
    return
end

function changeTau()
    return
end

function changeK()
    return
end

function changeExtTau()
    return
end

function changeExtK()
    return
end

export printStatus, changeOrder, changeTau, changeK, changeExtTau, changeExtK, eval, test

end
