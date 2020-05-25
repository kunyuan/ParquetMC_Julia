module Ver4Test
include("../parameter.jl")
include("propagator.jl")
include("vertex4.jl")
using .Vertex4
using .Propagator: interaction, green

function evalOneLoopVer4(chan, curr)
    weight = zero(VerWeight)

    ver4 = Vertex4.Ver4(0, 1, chan, 1, RIGHT, false, true)
    Vertex4.eval(ver4, curr.K[INL], curr.K[OUTL], curr.K[INR], curr.K[OUTR], 5, curr, true)
    weight = sum(ver4.weight)
    # for w in ver4.weight
    #     weight .+= w
    # end
    return weight
end

function testOneLoopVer4(curr)
    println("Testing ...")
    varK, varT = (curr.K, curr.T)
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
    weight[DI] = Lver[DI] * Rver[DI] * SPIN + Lver[EX] * Rver[DI] + Lver[DI] * Rver[EX]
    weight[EX] = Lver[EX] * Rver[EX]
    weight .*= gweight * PhaseFactor * Vertex4.SymFactor[T]

    cweight = zero(VerWeight)
    cweight[DI] = gweightbox * Lver[DI] * Rver[DI] * Lambda / (8.0 * pi) / Nf * SPIN
    cweight .*= PhaseFactor * Vertex4.SymFactor[T]
    testWeight = weight - cweight

    refweight = evalOneLoopVer4([T, TC], curr)
    @assert isapprox(testWeight, refweight, rtol = 1.0e-16) "\n$testWeight != $refweight\nGG: $gweight, GG_counter: $gweightbox\nLver: $Lver, Rver: $Rver"

    ### U and UC #####
    Ku = outR + K - inL
    Lver = interaction(inL, outR, Ku, K, false)
    Rver = interaction(K, Ku, inR, outL, false)

    dTau = varT[2] - varT[1]
    gweight = green(dTau, K) * green(-dTau, Ku)
    gweightbox = green(dTau, K) * green(-dTau, K)

    weight[DI] = Lver[EX] * Rver[EX]
    weight[EX] = Lver[DI] * Rver[DI] * SPIN + Lver[EX] * Rver[DI] + Lver[DI] * Rver[EX]
    weight .*= gweight * PhaseFactor * Vertex4.SymFactor[U]

    cweight = zero(VerWeight)
    cweight[EX] = gweightbox * Lver[DI] * Rver[DI] * Lambda / (8.0 * pi) / Nf * SPIN
    cweight .*= PhaseFactor * Vertex4.SymFactor[U]
    testWeight = weight - cweight

    refweight = evalOneLoopVer4([U, UC], curr)
    @assert isapprox(testWeight, refweight, rtol = 1.0e-16) "\n$testWeight != $refweight\nGG: $gweight, GG_counter: $gweightbox\nLver: $Lver, Rver: $Rver"

    Ks = inR + inL - K
    Lver = interaction(inL, Ks, inR, K, false)
    Rver = interaction(K, outL, Ks, outR, false)
    dTau = varT[2] - varT[1]
    gweight = green(dTau, K) * green(dTau, Ks)

    weight[DI] = Lver[EX] * Rver[DI] + Lver[DI] * Rver[EX]
    weight[EX] = Lver[DI] * Rver[DI] + Lver[EX] * Rver[EX]

    weight *= gweight * PhaseFactor * Vertex4.SymFactor[S]
    weight *= cos(2.0 * pi / Beta * dTau * 2.0)

    refweight = evalOneLoopVer4([S], curr)
    @assert isapprox(weight, refweight, rtol = 1.0e-16) "\n$weight != $refweight\nGG: $gweight, GG_counter: $gweightbox\nLver: $Lver, Rver: $Rver"
end

export testOneLoopVer4
end
