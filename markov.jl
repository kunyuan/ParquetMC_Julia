module Markov
include("parameter.jl")
include("grid.jl")
include("diag/vertex4.jl")
using Random
import .Vertex4
const UpdateNum = 9
const LastTidx = 2 * Order

function init(_counter, _rng)
    global rng = _rng
    global counter = _counter

    global currOrder = 0
    global currAbsWeight = 0.0
    global currScale, currExtTidx, currExtKidx, currExtAngidx = (1, 1, 1, 1)

    global varTau = rand(rng, LastTidx) * Beta
    varTau[LastTidx] = Grid.tau.grid[currExtTidx]
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

    global vertex4 = Vector{Vertex4.Ver4}(undef, 0)
    chan = [I, T, U, S, TC, UC]
    for order = 1:Order
        push!(vertex4, Vertex4.Ver4(0, order, chan, 5, 1, RIGHT, false))
    end
end

function printStatus()
    # Var.counter += 1
    ReWeight[2] = 1.6
    println(currExtAngidx)
    println(Beta)
end

function evaluate()
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

export printStatus, changeOrder, changeTau, changeK, changeExtTau, changeExtK, evaluate

end
