module Markov
include("parameter.jl")
include("grid.jl")
include("diag/vertex4.jl")
include("diag/polar.jl")
include("diag/vertex4_test.jl")
using Random, StaticArrays
import .Vertex4, .Ver4Test, .Polar
const UpdateNum = 9

# const curr = Main.State()
# const varT = Vector{Float}(undef, LastTidx)
# const varK = Vector{Mom}(undef, LastKidx)
const curr = Main.Curr
const rng = Main.Curr.rng
const varK = Main.Curr.K
const varT = Main.Curr.T

const ver4 = Vector{Vertex4.Ver4}(undef, 0)

# function init(_counter, _rng)
function init()

    #######  initialize MC variables  ################################

    ###### initialized diagram trees #######################################

    if DiagType == GAMMA
        chan = [I, T, U, S, TC, UC]
        for order = 1:Order
            push!(ver4, Vertex4.Ver4(0, order, chan, 1, RIGHT, false, true))
            GC.gc() # collect garabage, reduce memory usage
        end
        # Vertex4.visualize(vertex4[2])

    elseif DiagType == POLAR
        global polar = Vector{Polar.Polarization}(undef, 0)
        for o = 0:Order
            push!(polar, Polar.Polarization(o))
        end
    else
        throw("Not implemented!")
    end
end

function test()
    if DiagType == GAMMA
        Ver4Test.testOneLoopVer4(curr)
    end
end

function benchmark(o)
        # @allocated begin
    if DiagType == GAMMA
        @fastmath Vertex4.eval(
            ver4[o],
            varK[INL],
            varK[OUTL],
            varK[INR],
            varK[OUTR],
            5,
            true,
        )
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
