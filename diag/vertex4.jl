module Vertex4
include("../parameter.jl")
include("propagator.jl")
using .Propagator

function init(_varT::Vector{Int}, _varK::Vector{Mom})
    global varT = _varK
    global varK = _varK
end

struct Green
    Tidx::Vector{Tuple{Int,Int}}
    weight::Vector{Real}
    Green() = new([], [])
end

function eval(G::Green, K::Mom, IsAnomal = false)
    for (i, Tidx) in enumerate(G.Tidx)
        G.weight[i] = green(varT[Tidx[OUT]] - varT[Tidx[IN]], K)
    end
end

function addTidx(_obj, _Tidx)
    for (i, Tidx) in enumerate(_obj.Tidx)
        if Tidx == _Tidx
            return i
        end
    end
    push!(_obj.Tidx, _Tidx)
    push!(_obj.weight, 0.0)
    return length(_obj.weight)
end

struct Ver4
    level::Int
    loopNum::Int
    chan::Vector{Int} # list of channels
    Kidx::Int # intial Kidx
    TinLidx::Int # inital Tidx
    side::Int # right side vertex is always a full gamma4
    inBox::Bool

    G::Vector{Green}
    Tidx::Vector{Tuple{Int,Int,Int,Int}}
    weight::Vector{Real}

    function Ver4(lvl, loopNum, chan, kidx, tidx, side, inbox)
        g = @SVector [Green() for i = 1:6]
        ver4 = new(lvl, loopNum, chan, kidx, tidx, side, inbox, g, [], [])
        if loopNum == 0
            addTidx(ver4, (tidx, tidx, tidx, tidx))
        else
            UST = [c for c in chan if c != I]
            II = [c for c in chan if c == I]
            println(UST)
            println(II)
            for c in UST
                for ol = 0:loopNum - 1
                    bubble = Bubble(ver4, c, ol)
                end
            end
        end
        return ver4
    end
end

struct IdxMap
    Lver::Int
    Rver::Int
    G0::Int
    G::Int
    ver::Int
end

struct Bubble
    chan::Int
    Lver::Ver4
    Rver::Ver4
    # map::Vector{IdxMap}

    function Bubble(ver4::Ver4, chan::Int, oL::Int)
        FULL = [I, T, U, S, TC, UC]
        F = [I, U, S, TC, UC]
        V = [I, T, U, TC, UC]
        FULL_CT = [I, T, TC]
        F_CT = [I, TC]

        @assert chan != I "buildBubble can not process I channel!"
        @assert oL < ver4.loopNum "LVer loopNum must be smaller than the ver4 loopNum"

        lvl = ver4.level + 1
        oR = ver4.loopNum - 1 - oL
        TinLidx = ver4.TinLidx
        TinRidx = ver4.TinLidx + (oL + 1)
        KLidx = ver4.loopNum + 1
        KRidx = ver4.loopNum + 1 + oL

        if chan == T
            if ver4.inBox == false
                global Lver = Ver4(lvl, oL, F, KLidx, TinLidx, LEFT, false)
                global Rver = Ver4(lvl, oR, FULL, KRidx, TinRidx, RIGHT, false)
            else
                global Lver = Ver4(lvl, oL, F_CT, KLidx, TinLidx, LEFT, true)
                global Rver = Ver4(lvl, oR, FULL_CT, KRidx, TinRidx, RIGHT, true)
            end
        end

        return new(chan, Lver, Rver)

    end
end

    export init, eval, Ver4

end
