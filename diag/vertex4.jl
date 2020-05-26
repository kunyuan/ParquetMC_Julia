module Vertex4
include("../parameter.jl")
include("propagator.jl")

using StaticArrays

const ChanMap = [I, T, U, S, T, U]
const SymFactor = [1.0, -1.0, 1.0, -0.5, 1.0, -1.0]
const varT = Main.Curr.T
const varK = Main.Curr.K

const Span SVector{2, Int} #head, tail

mutable struct Green
    K::Int
    Tpair::SVector{2, Int}
    weight::double
end

struct Bubble
    loopNum::Int
    chan::Int
    Tidx::Int
    Kidx::Int
    inBox::Bool
    legK::SVector{4,Int}

    mapidx::Vector{Int} # map idx head and tail
end

struct Map
    Lver::Int
    Rver::Int
    G0idx::Int # shared G idx
    Gcidx::Int # channel dependent G idx
    idx::Int  # weight idx
end

struct VerTree
    bubble::Vector{Bubble}
    map::Vector{Map}(undef, 0)
    weight::Vector{VerWeight}(undef, 0) # weight of each ver4 in VerPool
    Tpair::Vector{SVector{4,Int}}(undef, 0)
    Ver4Pool()=new([], [], [])
end

const KPool =@Mector [Vector{Mom}(undef, 0) for o in 1:Order]
const GPool=@Mector [Vector{Green}(undef, 0) for o in 1:Order]
const VerPool = @MVector [Vector{Ver4Pool}() for o in 1:Order]

isequal(g1::Green, g2::Green) = (g1.K == g2.K&&g1.Tpair==g2.Tpair)
isequal(v1::Ver4, v2::Ver4) = (v1.loopNum == v2.loopNum &&v1.inBox==v2.inBox&& v1.Kidx == v2.Kidx && v1.chan == v2.chan&&v1.Tidx==v2.Tidx)
function isequal(k1::Mom, k2::Mom)
    # test isequal with two critera, just in case
    flag1 = prod([isapprox(k1[i], k2[i], rtol = 1.0e-10) for i in 1:DIM])
    flag2 = prod([isapprox(k1[i], k2[i], rtol = 1.0e-12) for i in 1:DIM])
    @assert flag1 == flag2 "Not sure if k1==k2!"
    return flag1
end

function pushObj!(Pool, Obj)
    # println(Obj)
    for (idx, _obj) in enumerate(Pool)
        _obj == Obj && return idx
    end
    push!(Pool, Obj)
    println("Pool $(length(Pool)), add: $Obj")
    return length(Pool)
end

function SetKchan!(legK, K, chan::Int, gK)
    KinL, KoutL, KinR, KoutR = legK
    if chan == T || chan == TC
        gK .= KoutL .+ K .- KinL
    elseif chan == U || chan == UC
        # can not be in box!
        gK .= KoutR .+ K .- KinL
    else
        # S channel, and cann't be in box!
        gK .= KinL .+ KinR .- K
    end
end

function GetLegK(legK, K0, Kc, chan::Int)
    KinL, KoutL, KinR, KoutR = legK
    if chan == T || chan == TC
        LLegK = (KinL, KoutL, Kc, K0)
        RLegK = (K0, Kc, KinR, KoutR)
    elseif chan == U || chan == UC
        LLegK = (KinL, KoutR, Kc, K0)
        RLegK = (K0, Kc, KinR, KoutL)
    else
        # S channel
        LLegK = (KinL, Kc, KinR, K0)
        RLegK = (K0, KoutL, Kc, KoutR)
    end
    return LLegK, RLegK
end

function ver(pool, level, loopNum, chan, legK, loopKidx, tidx, side, inbox)
    @assert length(legK) == 4 "4 leg K are expected!"
    # first push four legs into the 
    legKidx = @SVector [pushObj!(KPool, _k) for _k in legK]

    if loopNum <= 0
        # negative loopNum should never be used in evaluation
        Tpairidx = pushObj!(pool, (tidx, tidx, tidx, tidx))
        veridx = pushObj!(0, 0, 0, 0, inBox, legKidx, [])
        return ver4idx
    else
        UST = [c for c in chan if c != I]
        for c in UST
            for ol = 0:loopNum - 1
                map=bubble(pool, level, loopNum, ol, c, legK, loopKidx, tidx, side, inbox)
                veridx=pushObj!(loopNum, c, tidx, loopKidx, inbox, legKidx, map)
            end
        end
    end
end

function bubble(pool, level, loopNum, ol, chan, legK, loopKidx, tidx, side, inbox)
    @assert chan != I "bubble can not be I channel"
    @assert length(legK) == 4 "4 leg K are expected!"
    K0 = varK[loopKidx]
    K0idx = push!(KPool, K0) # the common K shared by all channel

    Kc = zero(Mom) # K of the corresponding channel
    SetKchan!(legK, K, c, Kc)
    Kcidx = push!(KPool, Kc)

    LLegK, RLegK = GetLegK(legK, K0)

    FULL = [I, T, U, S, TC, UC]
    F = [I, U, S, TC, UC]
    V = [I, T, U, TC, UC]
    FULL_CT = [I, T, TC]
    F_CT = [I, TC]

    lvl = level + 1
    oR = loopNum - 1 - oL
    TinLidx = Tidx
    TinRidx = Tidx + (oL + 1)
    Llopidx = loopKidx + 1
    Rlopidx = loopKidx + 1 + ol;

# function verTree(level, loopNum, chan, legK, loopKidx, tidx, side, inbox)
    if chan == T
        if ver4.inBox == false
            Lver = ver(lvl, oL, F, LLegK, Llopidx, TinLidx, LEFT, false)
            Rver = ver(lvl, oR, FULL, RLegK, Rlopidx, TinRidx, RIGHT, false)
        else
            Lver = ver(lvl, oL, F_CT, LLegK, Llopidx, TinLidx, LEFT, true)
            Rver = ver(lvl, oR, FULL_CT, RLegK, Rlopidx, TinRidx, RIGHT, true)
        end
    elseif chan == U
        @assert ver4.inBox == false "Ver4 in box can't have U diagrams!"
        Lver = ver(lvl, oL, F, LLegK, Llopidx, TinLidx, LEFT, false)
        Rver = vere(lvl, oR, FULL, RLegK, Rlopidx, TinRidx, RIGHT, false)
    elseif chan == S
        @assert ver4.inBox == false "Ver4 in box can't have S diagrams!"
        Lver = ver(lvl, oL, V, LLegK, Llopidx, TinLidx, LEFT, false)
        Rver = ver(lvl, oR, FULL, RLegK, Rlopidx, TinRidx, RIGHT, false)
    elseif chan == TC || chan == UC
        Lver = ver(lvl, oL, F_CT, LLegK, Llopidx, TinLidx, LEFT, true)
        Rver = ver(lvl, oR, FULL_CT, RLegK, Rlopidx, TinRidx, RIGHT, true)
    end

    # ############## construct IdxMap ########################################
    # map = []
    # G = ver4.G
    # for (lt, LvT) in enumerate(Lver.Tpair)
    #     for (rt, RvT) in enumerate(Rver.Tpair)
    #         GT0idx = addTidx(G[1], (LvT[OUTR], RvT[INL]))
    #         GTxidx, VerTidx = (0, 0)

    #         if chan == T
    #             VerTidx = addTidx(ver4, (LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR]))
    #             GTxidx = addTidx(G[T], (RvT[OUTL], LvT[INR]))
    #         elseif chan == U
    #             VerTidx = addTidx(ver4, (LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL]))
    #             GTxidx = addTidx(G[U], (RvT[OUTL], LvT[INR]))
    #         elseif chan == S
    #             VerTidx = addTidx(ver4, (LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR]))
    #             GTxidx = addTidx(G[S], (LvT[OUTL], RvT[INR]))
    #         elseif chan == TC || chan == UC
    #             # counterterms are equal-time
    #             VerTidx = addTidx(ver4, (LvT[INL], LvT[INL], LvT[INL], LvT[INL]))
    #             GTxidx = -1 # do not need it
    #         # GTxidx = addTidx(ver4.G[chan], (RvT[OUTL], LvT[INR]))
    #         else
    #             throw("This channel is invalid!")
    #         end

    #         @assert ver4.Tpair[end][1] == ver4.Tidx "InL Tidx must be the same for all Tpairs in the vertex4"
    #         @assert GT0idx != 0 && VerTidx != 0 "index cann't be zero!"
    #         push!(map, IdxMap(GT0idx, GTxidx, VerTidx))
    #     end
    # end
    # return new(chan, Lver, Rver, map)
    return 0, 0
end

end
