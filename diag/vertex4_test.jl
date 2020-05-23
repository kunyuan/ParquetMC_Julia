module Ver4Test
include("../parameter.jl")
include("propagator.jl")
include("vertex4.jl")
using .Vertex4
using .Propagator: interaction, green

# function init(_varT::Vector{Float}, _varK::Vector{Mom})
#     global varT = _varT
#     global varK = _varK
#     Vertex4.init(varT, varK)
# end

function evalOneLoopVer4(chan, varT, varK)
    weight = zero(VerWeight)

    ver4 = Vertex4.Ver4(0, 1, chan, 1, RIGHT, false, true)
    Vertex4.eval(ver4, varK[INL], varK[OUTL], varK[INR], varK[OUTR], 5, varT, varK, true)
    weight = sum(ver4.weight)
    # for w in ver4.weight
    #     weight .+= w
    # end
    return weight
end

function testOneLoopVer4(varT, varK)
    println("Testing ...")
    inL, outL, inR, outR = (varK[INL], varK[OUTL], varK[INR], varK[OUTR])
    K = varK[5]

    # T and TC
    Kt = outL + K - inL
    println(K)
    println(Kt)
    println(inR)
    println(outR)
    Lver = interaction(inL, outL, Kt, K, false)
    Rver = interaction(K, Kt, inR, outR, false)
    dTau = varT[2] - varT[1]
    gweight = green(dTau, K) * green(-dTau, Kt)
    gweightbox = green(dTau, K) * green(-dTau, K)

    weight = zero(VerWeight)
    weight.dir =
        Lver.dir * Rver.dir * SPIN + Lver.ex * Rver.dir + Lver.dir * Rver.ex
    weight.ex = Lver.ex * Rver.ex
    weight .*= gweight * PhaseFactor * Vertex4.SymFactor[T]

    cweight = zero(VerWeight)
    cweight.dir = gweightbox * Lver.dir * Rver.dir * Lambda / (8.0 * pi) / Nf * SPIN
    cweight .*= PhaseFactor * Vertex4.SymFactor[T]
    testWeight = weight - cweight

    refweight = evalOneLoopVer4([T, TC], varT, varK)

    if abs(testWeight.dir - refweight.dir) > 1.0e-16 || 
       abs(testWeight.ex - refweight.ex) > 1.0e-16
        printstyled(
                "$testWeight != $refweight\n",
                color = :red,
            )
        println("GG: $gweight, GG_counter: $gweightbox")
        println("Lver: $Lver")
        println("Rver: $Rver")
    end

    ### U and UC #####
    Ku = outR + K - inL
    Lver = interaction(inL, outR, Ku, K, false)
    Rver = interaction(K, Ku, inR, outL, false)

    dTau = varT[2] - varT[1]
    gweight = green(dTau, K) * green(-dTau, Ku)
    gweightbox = green(dTau, K) * green(-dTau, K)

    weight.dir = Lver.ex * Rver.ex
    weight.ex =
        Lver.dir * Rver.dir * SPIN + Lver.ex * Rver.dir + Lver.dir * Rver.ex
    weight .*= gweight * PhaseFactor * Vertex4.SymFactor[U]

    cweight = zero(VerWeight)
    cweight.ex = gweightbox * Lver.dir * Rver.dir * Lambda / (8.0 * pi) / Nf * SPIN
    cweight .*= PhaseFactor * Vertex4.SymFactor[U]
    testWeight = weight - cweight

    refweight = evalOneLoopVer4([U, UC], varT, varK)

    if abs(testWeight.dir - refweight.dir) > 1.0e-16 || 
       abs(testWeight.ex - refweight.ex) > 1.0e-16
        printstyled(
                    "$testWeight != $refweight\n",
                    color = :red,
                )
        println("GG: $gweight, GG_counter: $gweightbox")
        println("Lver: $Lver")
        println("Rver: $Rver")
    end

    Ks = inR + inL - K
    Lver = interaction(inL, Ks, inR, K, false)
    Rver = interaction(K, outL, Ks, outR, false)
    dTau = varT[2] - varT[1]
    gweight = green(dTau, K) * green(dTau, Ks)

    weight.dir = Lver.ex * Rver.dir + Lver.dir * Rver.ex
    weight.ex = Lver.dir * Rver.dir + Lver.ex * Rver.ex

    weight *= gweight * PhaseFactor * Vertex4.SymFactor[S]
    weight *= cos(2.0 * pi / Beta * dTau * 2.0)

    refweight = evalOneLoopVer4([S,], varT, varK)

    if abs(weight.dir - refweight.dir) > 1.0e-16 || 
       abs(weight.ex - refweight.ex) > 1.0e-16
        printstyled(
                    "$weight != $refweight\n",
                    color = :red,
                )
        println("GG: $gweight")
        println("Lver: $Lver")
        println("Rver: $Rver")
    end
end

export testOneLoopVer4
end
