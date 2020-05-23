module Ver4Test
include("../parameter.jl")
include("propagator.jl")
include("vertex4.jl")
using .Vertex4
using .Propagator: interaction, green

function init(_varT::Vector{Float}, _varK::Vector{Mom})
    global varT = _varT
    global varK = _varK
    Vertex4.init(varT, varK)
end

function evalOneLoopVer4(chan)
    weight = zero(VerWeight)

    ver4 = Vertex4.Ver4(0, 1, chan, 1, RIGHT, false, true)
    Vertex4.eval(ver4, varK[INL], varK[OUTL], varK[INR], varK[OUTR], 5, true)
    weight = sum(ver4.weight)
    # for w in ver4.weight
    #     weight .+= w
    # end
    return weight
end

function testOneLoopVer4()
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
    weight[DIR] =
        Lver[DIR] * Rver[DIR] * SPIN + Lver[EX] * Rver[DIR] + Lver[DIR] * Rver[EX]
    weight[EX] = Lver[EX] * Rver[EX]
    weight .*= gweight * PhaseFactor * Vertex4.SymFactor[T]

    cweight = zero(VerWeight)
    cweight[DIR] = gweightbox * Lver[DIR] * Rver[DIR] * Lambda / (8.0 * pi) / Nf * SPIN
    cweight .*= PhaseFactor * Vertex4.SymFactor[T]
    testWeight = weight - cweight

    refweight = evalOneLoopVer4([T, TC])

    if abs(testWeight[DIR] - refweight[DIR]) > 1.0e-10 || 
       abs(testWeight[EX] - refweight[EX]) > 1.0e-10
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

    weight[DIR] = Lver[EX] * Rver[EX]
    weight[EX] =
        Lver[DIR] * Rver[DIR] * SPIN + Lver[EX] * Rver[DIR] + Lver[DIR] * Rver[EX]
    weight .*= gweight * PhaseFactor * Vertex4.SymFactor[U]

    cweight = zero(VerWeight)
    cweight[EX] = gweightbox * Lver[DIR] * Rver[DIR] * Lambda / (8.0 * pi) / Nf * SPIN
    cweight .*= PhaseFactor * Vertex4.SymFactor[U]
    testWeight = weight - cweight

    refweight = evalOneLoopVer4([U, UC])

    if abs(testWeight[DIR] - refweight[DIR]) > 1.0e-10 || 
       abs(testWeight[EX] - refweight[EX]) > 1.0e-10
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

    weight[DIR] = Lver[EX] * Rver[DIR] + Lver[DIR] * Rver[EX]
    weight[EX] = Lver[DIR] * Rver[DIR] + Lver[EX] * Rver[EX]

    weight *= gweight * PhaseFactor * Vertex4.SymFactor[S]
    weight *= cos(2.0 * pi / Beta * dTau * 2.0)

    refweight = evalOneLoopVer4([S,])

    if abs(weight[DIR] - refweight[DIR]) > 1.0e-10 || 
       abs(weight[EX] - refweight[EX]) > 1.0e-10
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
