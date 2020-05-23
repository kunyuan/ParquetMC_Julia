module Polar
include("../parameter.jl")
include("vertex4.jl")
include("propagator.jl")

using .Vertex4: Green, Ver4
using .Propagator: green

function init(_varT::Vector{Float}, _varK::Vector{Mom})
    global varT = _varT
    global varK = _varK
end

struct Polarization
    order::Int
    ver4::Ver4
    G::Vector{Green}
    Gidx::Vector{NTuple{4,Int}}

    function Polarization(order::Int)
        # order<=1 should never be used in actual evaluation
        @assert InterTauNum(order) + 1 < LastTidx "More Tau variables are needed!"

        ver4 = Ver4(
            0,
            order - 2, # loopNum
            [I, T, U, S, TC, UC],
            2, # InL Tidx
            RIGHT,
            false, # no fast eval is needed
        )

        polar = new(order, ver4, [Green() for i = 1:4], [])
        for t in ver4.Tpair
            inL = Vertex4.addTidx(G[INL], (1, t[INL]))
            outL = Vertex4.addTidx(G[OUTL], (t[OUTL], 1))
            inR = Vertex4.addTidx(G[INR], (LastTidx, t[INR]))
            outR = Vertex4.addTidx(G[OUTR], (t[OUTR], LastTidx))
            push!(polar.Gidx, (inL, outL, inR, outR))
        end

        for g in polar.G
            for t in g.Tpair
                push!(g.weight, 0.0)
            end
        end

        return polar

    end
end

function eval(polar::Polarization)
    if polar.order == 0
        return 1.0
    elseif polar.order == 1
        Tau = varT[LastTidx] - varT[1]
        gWeight = green(Tau, varK[2]) * green(-Tau, varK[2] - varK[1])
        return -SPIN * gWeight * PhaseFactor
    end

    # order>=2
    K = [varK[2], varK[2] - varK[1], varK[3], varK[3] - varK[1]]
    G = polar.G
    Vertex4.eval.(G, K)
    Vertex4.eval(polar.ver4, K[INL], K[OUTL], K[INR], K[OUTR], 4)

    w = 0.0
    weight = polar.ver4.weight
    for (i, gidx) in polar.Gidx
        temp = (SPIN^2 * weight[i].dir + SPIN * weight[i].ex)
        for gidx in polar.Gidx
            temp *= G.weight[gidx]
        end
        w += temp
    end
    return w * PhaseFactor^2
end

export Polar, eval
end
