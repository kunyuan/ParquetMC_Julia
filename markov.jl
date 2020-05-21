module Markov
include("parameter.jl")
include("grid.jl")
using Random
const UpdateNum = 9

const LastTidx = 2 * Order

currOrder = 0
currAbsWeight = 0.0
currExtTidx = 1
currExtKidx = 1
currExtAngidx = 1

function init(_counter, _rng)
    global rng = _rng
    global counter = _counter
    global varTau = rand(rng, LastTidx) * Beta
    varTau[LastTidx] = Grid.tau.grid[currExtTidx]
    global varK = Vector{Mom}(undef, Order + 1)
    rand!(rng, varK)
    varK *= Kf
    println(Grid.tau.size)
end

function printStatus()
    # Var.counter += 1
    ReWeight[2] = 1.6
    println(Beta)
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

function evaluate()
    return
end

export printStatus, changeOrder, changeTau, changeK, changeExtTau, changeExtK, evaluate

end
