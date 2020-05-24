module Markov
include("parameter.jl")
include("grid.jl")
include("diag/vertex4.jl")
include("diag/polar.jl")
include("diag/vertex4_test.jl")
using Random, StaticArrays
import .Vertex4, .Ver4Test, .Polar
const UpdateNum = 9

const curr = State()
# const varT = Vector{Float}(undef, LastTidx)
# const varK = Vector{Mom}(undef, LastKidx)

const ver4 = Vector{Vertex4.Ver4}(undef, 0)

# function init(_counter, _rng)
function init()

    # global rng = _rng
    # global counter = _counter

    rng = MersenneTwister()
    counter = Counter()

    #######  initialize MC variables  ################################

    curr.T .= rand(rng, LastTidx) * Beta
    curr.T[LastTidx] = Grid.tau.grid[curr.extTidx]
    curr.T[1] = 0.0
    curr.T[2] = Beta / 2.0
    curr.T[3] = Beta / 3.0
    curr.T[4] = Beta / 4.0
    curr.T[5] = Beta / 5.0
    curr.T[6] = Beta / 6.0

    rand!(rng, curr.K)
    curr.K .*= Kf
    k0 = curr.K[1]

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
        k = [Grid.K.grid[curr.extKidx], 0.0, 0.0]
        curr.K[1] = k[1:DIM]
    end

    ###### initialized diagram trees #######################################
    # Vertex4.init(varT, varK)
    # Polar.init(varT, varK)
    # Ver4Test.init(varT, varK)

    if DiagType == GAMMA
        chan = [I, T, U, S, TC, UC]
        for order = 1:Order
            push!(ver4, Vertex4.Ver4(0, order, chan, 1, RIGHT, false, true))
            GC.gc() # collect garabage, reduce memory usage
        end
        # Vertex4.visualize(vertex4[2])

        for o = 1:Order
            println("Order $o")
            for i = 1:10
                Vertex4.eval(
                    ver4[o],
                    curr.K[INL],
                    curr.K[OUTL],
                    curr.K[INR],
                    curr.K[OUTR],
                    5,
                    curr,
                    true,
                )
                # println(sum(ver4[o].weight))
                # println((ver4[o].weight))
            end
        end

        for o = 1:Order
            @time begin
                # @allocated begin
                println("Order $o")
                for i = 1:1000
                    Vertex4.eval(
                        ver4[o],
                        curr.K[INL],
                        curr.K[OUTL],
                        curr.K[INR],
                        curr.K[OUTR],
                        5,
                        curr,
                        true,
                    )
                end
            end
        end
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
