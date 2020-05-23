module Markov
include("parameter.jl")
include("grid.jl")
include("diag/vertex4.jl")
include("diag/polar.jl")
include("diag/vertex4_test.jl")
using Random
import .Vertex4, .Ver4Test, .Polar
const UpdateNum = 9

function init(_counter, _rng)
    global rng = _rng
    global counter = _counter

    #######  initialize MC variables  ################################
    global currOrder = 0
    global currAbsWeight = 0.0
    global currScale, currExtTidx, currExtKidx, currExtAngidx = (1, 1, 1, 1)

    global varT = rand(rng, LastTidx) * Beta
    varT[LastTidx] = Grid.tau.grid[currExtTidx]

    global varK = Vector{Mom}(undef, Order + 4)
    rand!(rng, varK)
    varK *= Kf

    if DiagType == GAMMA
        kL = [Kf, 0.0, 0.0]
        θ = acos(Grid.angle.grid[currExtAngidx])
        kR = [Kf * cos(θ), Kf * sin(θ), 0.0]
        varK[INL], varK[OUTL] = (kL[1:DIM], kL[1:DIM], kR[1:DIM], kR[1:DIM])
    else
        k = [Grid.K.grid[currExtKidx], 0.0, 0.0]
        varK[1] = k[1:DIM]
    end

    ###### initialized diagram trees #######################################
    Vertex4.init(varT, varK)
    Polar.init(varT, varK)
    Ver4Test.init(varT, varK)

    if DiagType == GAMMA
        global vertex4 = Vector{Vertex4.Ver4}(undef, 0)
        chan = [I, T, U, S, TC, UC]
        for order = 1:Order
            push!(vertex4, Vertex4.Ver4(0, order, chan, 1, RIGHT, false, true))
        end
        # Vertex4.visualize(vertex4[end])
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
        Ver4Test.testOneLoopVer4()
    end
end

function printStatus()
    # Var.counter += 1
    ReWeight[2] = 1.6
    println(currExtAngidx)
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
