module Markov
include("parameter.jl")
include("grid.jl")
include("diag/vertex4.jl")
include("diag/polar.jl")
include("diag/vertex4_test.jl")
using Random, StaticArrays, Printf
import .Vertex4, .Ver4Test, .Polar
const UpdateNum = 6
const INCREASE_ORDER, DECREASE_ORDER, CHANGE_EXTTAU, CHANGE_EXTK, CHANGE_TAU, CHANGE_K =
    1:UpdateNum
const Name = [
    "increase_order",
    "decrease_order",
    "change_ExtTau",
    "change_ExtK",
    "change_Tau",
    "change_K",
]
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
    if rand(rng) < 0.5
        curr.order == Order && return # already at the highest order
        name = INCREASE_ORDER
        newOrder = curr.order + 1
        if curr.order == 0
            curr.extTidx, propT = createExtIdx(TauGridSize)
            varT[LastTidx] = Grid.tau.grid[curr.extTidx]
            curr.extKidx, propK = createExtIdx(KGridSize)
            varK[0] = Grid.K.grid(curr.extKidx)
        else
            varT[lastInnerTidx(curr.order)]
        end

    else

    end
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

@inline createExtIdx(GridSize) = rand(rng, 1:GridSize), GridSize
@inline removeExtIdx(GridSize) = 1.0 / GridSize
@inline shiftExtIdx(GridSize) = rand(rng, 1:GridSize), 1.0

# newTau, Prop
@inline createTau() = rand(rng) * Beta, Beta
@inline removeTau(oldTau) = 1.0 / Beta
@inline function shiftTau(oldTau)
    x = rand(rng)
    newTau = 0.0
    prop = 1.0
    if x < 1.0 / 3
        newTau = oldTau + 4 * Ef * (rand(rng) - 0.5), prop
    elseif x < 2.0 / 3
        newTau = -oldTau
    else
        newTau = rand(rng) * Beta
    end

    if newTau < 0.0
        newTau += Beta
    elseif newTau > Beta
        newTau -= Beta
    elseif newTau == 0.0 || newTau == 0.0
        prop = 0.0
    end
    return newTau, prop
end

@inline function createK(oldK)
    dK = Kf / 2.0
    Kamp = Kf + (rand(rng) - 0.5) * 2.0 * dK
    Kamp <= 0.0 && return oldK, 0.0
    newK = Mom()
    prop = 0.0
    # Kf-dK<Kamp<Kf+dK 
    ϕ = 2.0 * pi * rand(rng)
    if DIM == 3
        θ = pi * rand(rng)
        newK[1] = Kamp * cos(ϕ) * sin(θ)
        newK[2] = Kamp * sin(ϕ) * sin(θ)
        newK[3] = Kamp * cos(θ)
        prop = (2.0 * dK) * (2.0 * pi) * pi * (sin(θ) * Kamp^2)
        # prop density of KAmp in [Kf-dK, Kf+dK), prop density of Phi
        # prop density of Theta, Jacobian
    else
        # DIM==2
        newK[1] = Kamp * cos(θ)
        newK[2] = Kamp * sin(θ)
        prop = (2.0 * dK) * (2.0 * pi) * (Kamp)
        # prop density of KAmp in [Kf-dK, Kf+dK), prop density of Phi, Jacobian
    end
    return newK, prop
end

@inline function removeK(oldK)
    dK = Kf / 2.0
    Kamp = norm(oldK)
    (Kamp < Kf - dK || Kamp > Kf + dK) && return 0.0
    if DIM == 3
        sinTheta = sqrt(oldK[1]^2 + oldK[2]^2) / Kamp
        if sinTheta > 1.0e-15
            return 1.0 / (2.0 * dK * 2.0 * pi^2 * sinTheta * Kamp^2)
        else
            return 0.0
        end
    else
        # DIM==2
        return 1.0 / (2.0 * dK * 2.0 * pi * Kamp)
    end
end

@inline function shiftK(oldK)
    x = rand(rng)
    if x < 1.0 / 3
        dK = Beta > 1.0 ? Kf / Beta * 3.0 : Kf
        return oldK + rand(rng, DIM) * dK, 1.0
    elseif x < 2.0 / 3
        λ = 1.5
        ratio = 1.0 / λ + rand(rng) * (λ - 1.0 / λ)
        prop = (DIM == 2) ? 1.0 : ratio
        return oldK * ratio, prop
    else
        return -oldK, 1.0
    end
end

const barbar = "====================================================================================="
const bar = "-------------------------------------------------------------------------------------"

function printStatus()
    # Var.counter += 1
    println(barbar)
    println("Step:", curr.step)
    println(bar)
    for i = 1:UpdateNum
        @printf("%-12s %12s %12s %12s\n", Name[i], "Proposed", "Accepted", "Ratio  ")
        for o = 1:Order
            @printf(
                "  Order%2d:   %12.6f %12.6f %12.6f\n",
                o,
                Proposed[i, o + 1],
                Accepted[i, o + 1],
                Accepted[i, o + 1] / Proposed[i, o + 1]
            )
        end
        println(bar)
    end
end

end
