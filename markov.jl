module Markov
include("parameter.jl")
include("grid.jl")
include("diag/vertex4.jl")
include("diag/polar.jl")
include("diag/vertex4_test.jl")
using Random, StaticArrays, Printf
import .Vertex4, .Ver4Test, .Polar
const UpdateNum = 6
const INCREASE_ORDER, DECREASE_ORDER, CHANGE_EXTTAU, CHANGE_EXTK, CHANGE_TAU, CHANGE_K = 1:UpdateNum
const Name = ["increase_order", "decrease_order", "change_ExtTau", "change_ExtK", "change_Tau", "change_K"]
const Accepted = (@MArray zeros(Float, UpdateNum, Order + 1)) .+ 1.0e-10
const Proposed = (@MArray zeros(Float, UpdateNum, Order + 1)) .+ 1.0e-10

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

function measure()
    return
end

function save()
    return
end

function reweight()
    return
end

function eval(order)
    if DiagType == POLAR
        return Polar[order].eval()
    end
end

function changeOrder()
    return
end

function changeTau()
    Proposed[CHANGE_TAU, curr.order + 1] += 1
    Accepted[CHANGE_TAU, curr.order + 1] += 1
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

const barbar = "====================================================================================="
const bar = "-------------------------------------------------------------------------------------"

function printStatus()
    # Var.counter += 1
    println(barbar)
    println("Step:", curr.step)
    println(bar)
    for i in 1:UpdateNum
        @printf("%-12s %12s %12s %12s\n", Name[i], "Proposed", "Accepted", "Ratio  ")
        for o in 1:Order
            @printf("  Order%2d:   %12.6f %12.6f %12.6f\n", o, Proposed[i, o + 1], Accepted[i, o + 1], Accepted[i, o + 1] / Proposed[i, o + 1])
        end
        println(bar)
    end
end

end
