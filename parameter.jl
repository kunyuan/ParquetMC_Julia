using StaticArrays
################ Global parameters  ##################################
const GAMMA, SIGMA, POLAR, DELTA = (1, 2, 3, 4)
const DiagType = SIGMA
# const DiagType = GAMMA
const DIM, SPIN = (3, 2)
const (BoldG, BoldVer4) = (true, true)
const Order = 2
const beta, Rs, Mass2, Lambda, maxK = (40.0, 1.0, 0.0, 1.0, 3.0)
const TauGridSize, KGridSize, AngGridSize = (128, 32, 32)
# MC parameter
const TotalStep = 101e6
const ReWeight = [1.0, 3.0, 30.0, 1.0, 0.2]
const PrintTime, SaveTime, ReWeightTime, MessageTime, CollectTime = (10, 10, 30, 10, 10)

############  Derived parameters ###################################
@assert DIM == 2 || DIM == 3 "DIM not implemented!"
@assert length(ReWeight) > Order + 1 "More ReWeight elements are needed!"
const Kf = (DIM == 3) ? (9.0 * pi / 4.0)^(1.0 / 3) / Rs : sqrt(2.0) / Rs
const Ef, Mu, Nf = (Kf^2, Kf^2, Kf / (4.0 * pi^2) * SPIN)
const Beta = beta / Ef # rescale the temperature
const MaxK = maxK * Kf

############ Global External Variable Grid #########################
########### Other constants  #####################################
const Mom = SVector{3,Float64}
const IN, OUT = (1, 2)
const INL, OUTL, INR, OUTR = (1, 2, 3, 4)
const DIR, EX = (1, 2)
const DOWN, UP = (1, 2)
const LEFT, RIGHT = (1, 2)
const I, T, U, S, TC, UC = (1, 2, 3, 4, 5, 6)

struct VerWeight
    dir::Real
    ex::Real
    function VerWeight()
        new(0.0, 0.0)
    end
end

########## Global function  #######################################
function Counter()
    _counter::Int128 = 0
    step() = (_counter += 1)
    state() = _counter
    () -> (step, state) # make step() and state() public
end

@inline function InterTauNum(order)
    if DiagType == SIGMA || DiagType == DELTA
        return order - 2
    elseif DiagType == POLAR
        return order - 1
    else
        # GAMMA
        return order
    end
end

squaredNorm(k) = DIM == 3 ? k[1]^2 + k[2]^2 + k[3]^2 : k[1]^2 + k[2]^2
norm(k) = DIM == 3 ? sqrt(k[1]^2 + k[2]^2 + k[3]^2) : sqrt(k[1]^2 + k[2]^2)
dot(k, q) = DIM == 3 ? k[1] * q[1] + k[2] * q[2] + k[3] * q[3] : k[1] * q[1] + k[2] * q[2]
