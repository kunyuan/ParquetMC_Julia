include("utility/grid.jl")
using StaticArrays
using Random
using .Grid
################ Global parameters  ##################################
const global GAMMA, SIGMA, POLAR, DELTA = (1, 2, 3, 4)
const global DIM, SPIN = (3, 2)
const global BoldG = true
const global DiagType = SIGMA
const global Order = 2
const global Beta, Rs, Mass2, Lambda = (40.0, 1.0, 0.0, 1.0)
# Grid parameter
const global MaxK = 3.0
const global TauGridSize, KGridSize, AngGridSize = (128, 32, 32)
# MC parameter
const global TotalStep = 101e6
const global ReWeight = [1.0, 3.0, 30.0, 1.0, 0.2]
const global PrintTime, SaveTime, ReWeightTime, MessageTime, CollectTime =
    (10, 10, 30, 10, 10)
const global RNG = MersenneTwister(); # system random number is used for seeding
global PID = 0
global Counter = Int128(0)

############  Derived parameters ###################################
@assert DIM == 2 || DIM == 3 "DIM not implemented!"
@assert length(ReWeight) > Order + 1 "More ReWeight elements are needed!"
const global Kf = (DIM == 3) ? (9.0 * pi / 4.0)^(1.0 / 3) / Rs : sqrt(2.0) / Rs 
const global Ef, Mu, Nf = (Kf^2, Kf^2, Kf / (4.0 * pi^2) * SPIN)
Beta /= Ef # rescale the temperature
MaxK *= Kf

########### Global MC variables  ##################################

############ Global External Variable Grid #########################
# cos(theta) grids
const global AngGrid = Grid.UniformGrid(-1.0, 1.0, AngGridSize)

# Tau Grid construction
lambda = Beta / Ef / TauGridSize / 3.0
c1 = Coeff([0.0, Beta / 2.0], [1.0, TauGridSize / 2 + 0.5], lambda, true)
c2 = Coeff([Beta / 2.0, Beta], [TauGridSize / 2 - 0.5, TauGridSize], lambda, false)
const global TauGrid = Grid.LogGrid([c1, c2], [1, Int(TauGridSize / 2) + 1,TauGridSize])
TauGrid.grid[1], TauGrid.grid[end] = 1.0e-8, Beta - 1.0e-8

# MomGrid construction
if DiagType == SIGMA
    # Fermionic Grid
    @assert MaxK > Kf "MaxK must larger than Kf!"
    kFIdx = floor(Int, KGridSize * log(Kf) / log(MaxK - Kf))
    lambda = sqrt(Ef * Beta) / kFIdx
    c1 = Coeff([0.0, Kf], [1.0, kFIdx + 1.0], lambda, false)
    c2 = Coeff([Kf, MaxK], [kFIdx, KGridSize], lambda, true)
    const global KGrid = Grid.LogGrid([c1, c2], [1, kFIdx, KGridSize])
    KGrid.grid[1] = 1.0e-6
else
    # Bosonic Grid
    @assert MaxK > 2.0 * Kf "MaxK must larger than 2Kf!"
    kFIdx = Int(KGridSize / 3)
    twokFIdx = Int(KGridSize / 3 * 2)
    lambda = sqrt(Ef * Beta) / kFIdx
    c1 = Coeff([0.0, Kf], [1.0, kFIdx + 1.0], lambda, true)
    c2 = Coeff([Kf, 2.0 * Kf], [kFIdx, twokFIdx + 1.0], lambda, false)
    c3 = Coeff([2.0 * Kf, MaxK], [twokFIdx, KGridSize], lambda, true)
    const global KGrid = Grid.LogGrid([c1, c2, c3], [1, kFIdx, twokFIdx, KGridSize])
    KGrid.grid[1] = 1.0e-6
end



########### Other constants  #####################################
const global Mom = SVector{3,Float64}
const global IN, OUT = (1, 2)
const global INL, OUTL, INR, OUTR = (1, 2, 3, 4)
const global DIR, EX = (1, 2)
const global DOWN, UP = (1, 2)






