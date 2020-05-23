module Ver4Test
include("../parameter.jl")
include("propagator.jl")
include("vertex4.jl")
using .Vertex4
using .Propagator:interaction, green

function init(_varT::Vector{Float}, _varK::Vector{Mom})
    global varT = _varT
    global varK = _varK
    Vertex4.init(varT, varK)
end

function evalOneLoopVer4(chan)
    weight = VerWeight()
    
    ver4 = Vertex4.Ver4(0, 1, chan, 1, RIGHT, false, true);
    Vertex4.eval(ver4, varK[INL], varK[OUTL], varK[INR], varK[OUTR], 5,  true)
    for w in ver4.weight
        weight.dir += w.dir
        weight.ex += w.ex
    end
    return weight
end

function testOneLoopVer4()
    println("Testing ...")
    inL, outL, inR, outR = (varK[INL], varK[OUTL], varK[INR], varK[OUTR])
    K = varK[5]

    # T and TC
    Kt = outL + K - inL
    Lver = interaction(inL, outL, Kt, K, false)
    Rver = interaction(K, Kt, inR, outR, false)
    dTau = varT[2] - varT[1]
    gweight = green(dTau, K) * green(-dTau, Kt)
    gweightbox = green(dTau, K) * green(-dTau, K)

    weight = VerWeight()
    weight.dir = gweight * (Lver.dir * Rver.dir * SPIN + Lver.ex * Rver.dir + Lver.dir * Rver.ex)
    weight.ex = gweight * (Lver.ex * Rver.ex)
    weight *= PhaseFactor * Vertex4.SymFactor[T]

    cweight = VerWeight()
    cweight.dir = gweightbox * Lver.dir * Rver.dir * Lambda / (8.0 * pi) / Nf * SPIN
    cweight *= PhaseFactor * Vertex4.SymFactor[T]
    testWeight = weight - cweight

    refweight = evalOneLoopVer4([T, TC])

    if abs(testWeight.dir - refweight.dir) < 1.0e-10 && abs(testWeight.ex - refweight.ex)
        printstyled("($(testWeight.dir), $(testWeight.ex)) != ($(weight.dir), $(weight.ex))")
        println("GG: $gweight, GG_counter: $gweightbox")
        println("Lver: ($(Lver.dir), $(Lver.ex))")
        println("Rver: ($(Rver.dir), $(Rver.ex))")
    end
end

export testOneLoopVer4
end