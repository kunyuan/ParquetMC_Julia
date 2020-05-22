module Propagator
include("../parameter.jl")

function green(tau::Real, k::Mom, gtype::Int = 0, scale::Real = 0.0)
    if gtype == 1 || tau == 0.0
        tau = -1.0e-12
    end
    s = sign(tau)
    if tau < 0.0
        tau += Beta
    end
    Ek = squaredNorm(k)
    x = Beta * (Ek - Mu) * 0.5
    y = 2.0 * tau / Beta - 1.0
    if -100.0 < x < 100.0
        return s * exp(-x * y) / (2.0 * cosh(x))
    elseif x > 100.0
        return s * exp(-x * (y + 1.0))
    else
        return s * exp(x * (1.0 - y))
    end
end

function interaction(mom::Mom, verOrder::Int = 0)
    q2 = squaredNorm(mom)
    weight = q2 > 1.0e-14 ? -8.0 * pi / (q2 + Mass2 + Lambda) : 0.0
    if verOrder > 0
        weight *= (weight * Lambda / 8.0 / pi)^verOrder
    end
    return weight
end

function interaction(kInL::Mom, kOutL::Mom, kInR::Mom, kOutR::Mom, Boxed::Bool, extQ::Real)
    weight = VerWeight()
    qDi2 = squaredNorm(kInL - kOutL)
    weight.dir = -8.0 * pi / (qDi2 + Mass2 + Lambda)
    if (DiagType == SIGMA && qDi2 < 1.0e-14) ||
       (DiagType == POLAR && abs(qDi2 - extQ^2) < 1.0e-14)
        weight.dir = 0.0
    end

    if Boxed == false
        qEx2 = squaredNorm(kInL - kOutR)
        weight.ex = 8.0 * pi / (qEx2 + Mass2 + Lambda)
        if (DiagType == SIGMA && qEx2 < 1.0e-14) ||
           (DiagType == POLAR && abs(qEx2 - extQ^2) < 1.0e-14)
            weight.ex = 0.0
        end
    else
        weight.ex = 0.0
    end

end

counterBubble(K::Mom) =
    Lambda / (8.0 * pi * Kf) * green(Beta / 2.0, K) * green(-Beta / 2.0, K)

angle(K1::Mom, K2::Mom) = dot(K1, K2) / norm(K1) / norm(K2)

export green, interaction, angle, counterBubble
end
