module Markov
include("parameter.jl")
include("grid.jl")
include("diag/vertex4.jl")
include("diag/polar.jl")
include("diag/vertex4_test.jl")
include("observable.jl")
include("utility/utility.jl")
using Random, StaticArrays, Printf, Dates
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
const polar = Vector{Polar.Polarization}(undef, 0)

# function init(_counter, _rng)
function init()
    #######  initialize MC variables  ################################
    global oneBody = Observable.OneBody()
    println(oneBody.norm)
    # println(typeof(oneBody))
    ###### initialized diagram trees #######################################

    if DiagType == GAMMA
        chan = [I, T, U, S, TC, UC]
        for order = 1:Order
            push!(ver4, Vertex4.Ver4(0, order, chan, 1, RIGHT, false, true))
            GC.gc() # collect garabage, reduce memory usage
        end
        # Vertex4.visualize(vertex4[2])

    elseif DiagType == POLAR
        for o = 1:Order
            push!(polar, Polar.Polarization(o))
        end
    else
        throw("Not implemented!")
    end

end

@fastmath function test()
    if DiagType == GAMMA
        Ver4Test.testOneLoopVer4(curr)
    end
end

@fastmath function eval(order)
    order == 0 && return 1.0
    if DiagType == POLAR
        return Polar.eval(polar[order])
    elseif DiagType == GAMMA
        Vertex4.eval(ver4[order], varK[INL], varK[OUTL], varK[INR], varK[OUTR], 5, true)
        chanWeight = sum(ver4[order].weight)
        return chanWeight[DI] + chanWeight[EX] / SPIN
    else
        throw("Not implemented!")
    end
end

function measure()
    factor = 1.0 / curr.absWeight / ReWeight[curr.order + 1]
    if DiagType == POLAR || DiagType == SIGMA || DiagType == DELTA
        Observable.measure(oneBody, eval(curr.order), factor)
    end
end

function save()
    if DiagType == POLAR || DiagType == SIGMA || DiagType == DELTA
        Observable.save(oneBody)
    end
end

function reweight()
    return
end

function increaseOrder()
    curr.order == Order && return # already at the highest order
    newOrder = curr.order + 1
    if curr.order == 0
        # create new external Tau and K
        newextTidx, propT = createExtIdx(TauGridSize)
        varT[LastTidx] = Grid.tau.grid[curr.extTidx]
        newextKidx, propK = createExtIdx(KGridSize)
        varK[1][1] = Grid.K.grid[curr.extKidx]
    else
        # create new internal Tau and K
        varT[lastInnerTidx(curr.order)], propT = createTau()
        varK[lastInnerKidx(curr.order)], propK = createK()
    end
    prop = propT * propK

    newAbsWeight = abs(eval(newOrder))
    R = prop * newAbsWeight * ReWeight[newOrder + 1] / curr.absWeight / ReWeight[curr.order + 1]
    Proposed[INCREASE_ORDER, curr.order + 1] += 1
    if rand(rng) < R
        Accepted[INCREASE_ORDER, curr.order + 1] += 1
        curr.order = newOrder
        curr.absWeight = newAbsWeight
        if curr.order == 0
            curr.extTidx = newextTidx
            curr.extKidx = newextKidx
        end
    end
end

function decreaseOrder()
    curr.order == 0 && return
    newOrder = curr.order - 1
    if newOrder == 0
        # remove external Tau and K
        propT = removeExtIdx(TauGridSize)
        propK = removeExtIdx(KGridSize)
    else
        # remove internal Tau and K
        propT = removeTau()
        propK = removeK(varK[lastInnerKidx(curr.order)])
    end
    prop = propT * propK
    newAbsWeight = abs(eval(newOrder))
    R = prop * newAbsWeight * ReWeight[newOrder + 1] / curr.absWeight / ReWeight[curr.order + 1]
    Proposed[DECREASE_ORDER, curr.order + 1] += 1
    if rand(rng) < R
        Accepted[DECREASE_ORDER, curr.order + 1] += 1
        curr.order = newOrder
        curr.absWeight = newAbsWeight
    end
end

function changeTau()
    # Proposed[CHANGE_TAU, curr.order + 1] += 1
    # Accepted[CHANGE_TAU, curr.order + 1] += 1
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
@inline removeTau() = 1.0 / Beta
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

@inline function createK()
    dK = Kf / 2.0
    newK = zero(Mom)
    Kamp = Kf + (rand(rng) - 0.5) * 2.0 * dK
    Kamp <= 0.0 && return newK, 0.0
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
    printstyled(Dates.now(), color = :yellow)
    println("\nStep:", curr.step)
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
    println(progressBar(round(curr.step / 1000_000, digits = 2), TotalBlock))
    println()
end

end
