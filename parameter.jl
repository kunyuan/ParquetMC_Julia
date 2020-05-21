include("utility/grid.jl")
using Random
using StaticArrays
using .Grid
################ Global parameters  ##################################
const GAMMA, SIGMA, POLAR, DELTA = (1, 2, 3, 4)
const DIM, SPIN = (3, 2)
const BoldG = true
const DiagType = SIGMA
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
# cos(theta) grids
const AngGrid = Grid.UniformGrid(-1.0, 1.0, AngGridSize)

let
    # Tau Grid construction
    lambda = Beta / Ef / TauGridSize / 3.0
    c1 = Coeff([0.0, Beta / 2.0], [1.0, TauGridSize / 2 + 0.5], lambda, true)
    c2 = Coeff([Beta / 2.0, Beta], [TauGridSize / 2 - 0.5, TauGridSize], lambda, false)
    const global TauGrid =
        Grid.LogGrid([c1, c2], [1, Int(TauGridSize / 2) + 1, TauGridSize])
    TauGrid.grid[1], TauGrid.grid[end] = (1.0e-8, Beta - 1.0e-8)

    # MomGrid construction
    if DiagType == SIGMA
        # Fermionic Grid
        @assert MaxK > Kf "MaxK must larger than Kf!"
        kFi = floor(Int, KGridSize * log(Kf) / log(MaxK - Kf))
        lambda = sqrt(Ef * Beta) / kFi
        c1 = Coeff([0.0, Kf], [1.0, kFi + 1.0], lambda, false)
        c2 = Coeff([Kf, MaxK], [kFi, KGridSize], lambda, true)
        const global KGrid = Grid.LogGrid([c1, c2], [1, kFi, KGridSize])
        KGrid.grid[1] = 1.0e-6
    else
        # Bosonic Grid
        @assert MaxK > 2.0 * Kf "MaxK must larger than 2Kf!"
        kFi, twokFi = (Int(KGridSize / 3), Int(KGridSize / 3 * 2))
        lambda = sqrt(Ef * Beta) / kFi
        c1 = Coeff([0.0, Kf], [1.0, kFi + 1.0], lambda, true)
        c2 = Coeff([Kf, 2.0 * Kf], [kFi, twokFi + 1.0], lambda, false)
        c3 = Coeff([2.0 * Kf, MaxK], [twokFi, KGridSize], lambda, true)
        const global KGrid = Grid.LogGrid([c1, c2, c3], [1, kFi, twokFi, KGridSize])
        KGrid.grid[1] = 1.0e-6
    end
end

########### Other constants  #####################################
const Mom = SVector{3,Float64}
const IN, OUT = (1, 2)
const INL, OUTL, INR, OUTR = (1, 2, 3, 4)
const DIR, EX = (1, 2)
const DOWN, UP = (1, 2)
