module Propagator
include("../parameter.jl")

@inline function green(tau::Float, k::Mom, scale::Float = 0.0)
    if tau == 0.0 
        tau = -1.0e-12
    end
    s = sign(tau)
    if tau < 0.0
        tau += Beta
    end
    # @assert tau < Beta "tau must be [0.0, Beta)"
    Ek = squaredNorm(k) - Mu
    x = Beta * Ek * 0.5
    y = 2.0 * tau / Beta - 1.0
    if -100.0 < x < 100.0
        return s * exp(-x * y) / (2.0 * cosh(x))
    elseif x >= 100.0
        return s * exp(-x * (y + 1.0))
    else # x<=-100.0
        return s * exp(x * (1.0 - y))
    end
end

function interaction(mom::Mom, verOrder::Int = 0)
    q2 = squaredNorm(mom)
    weight = q2 > 1.0e-14 ? -8π / (q2 + Mass2 + Lambda) : 0.0
    if verOrder > 0
        weight *= (weight * Lambda / 8π)^verOrder
    end
    return weight
end

function interaction(
    kInL::Mom,
    kOutL::Mom,
    kInR::Mom,
    kOutR::Mom,
    Boxed::Bool,
    extQ::Real = -1.0,
)
    weight = VerWeight(0.0, 0.0)
    qDi2 = squaredNorm(kInL - kOutL)
    weight[DI] = -8π / (qDi2 + Mass2 + Lambda)
    if (DiagType == SIGMA && qDi2 < 1.0e-14) ||
       (DiagType == POLAR && abs(qDi2 - extQ^2) < 1.0e-14)
        weight[DI] = 0.0
    end

    if Boxed == false
        qEx2 = squaredNorm(kInL - kOutR)
        weight[EX] = 8π / (qEx2 + Mass2 + Lambda)
        if (DiagType == SIGMA && qEx2 < 1.0e-14) ||
           (DiagType == POLAR && abs(qEx2 - extQ^2) < 1.0e-14)
            weight.ex = 0.0
        end
    else
        weight[EX] = 0.0
    end
    return weight
end

counterBubble(K::Mom) =
    Lambda / (8π * Nf) * green(Beta / 2.0, K) * green(-Beta / 2.0, K)

angle(K1::Mom, K2::Mom) = dot(K1, K2) / norm(K1) / norm(K2)

end
