module Vertex4
include("../parameter.jl")
include("propagator.jl")
using .Propagator: interaction, green, counterBubble
using D3Trees

const ChanMap = [I, T, U, S, T, U]
const SymFactor = [1.0, -1.0, 1.0, -0.5, 1.0, -1.0]

function init(_varT::Vector{Float}, _varK::Vector{Mom})
    global varT = _varT
    global varK = _varK
end

struct Green
    K::Mom
    Tpair::Vector{Tuple{Int,Int}}
    weight::Vector{Float}
    Green() = new(zero(Mom), [], [])
end

function eval(G::Green, K::Mom, IsAnomal = false)
    G.K .= K
    for (i, t) in enumerate(G.Tpair)
        G.weight[i] = green(varT[t[OUT]] - varT[t[IN]], K)
    end
end

function addTidx(_obj, _Tidx)
    for (i, Tidx) in enumerate(_obj.Tpair)
        if Tidx == _Tidx
            return i
        end
    end
    push!(_obj.Tpair, _Tidx)
    return length(_obj.Tpair)
end

struct IdxMap
    Lver::Int
    Rver::Int
    G::Int
    Gx::Int
    ver::Int
end

struct Bubble{_Ver4} # template Bubble to avoid mutually recursive struct
    chan::Int
    Lver::_Ver4
    Rver::_Ver4
    map::Vector{IdxMap}

    function Bubble{_Ver4}(ver4::_Ver4, chan::Int, oL::Int) where {_Ver4}
        # printstyled("Bubble creation!\n", color = :red)
        FULL = [I, T, U, S, TC, UC]
        F = [I, U, S, TC, UC]
        V = [I, T, U, TC, UC]
        FULL_CT = [I, T, TC]
        F_CT = [I, TC]

        @assert chan != I "buildBubble can not process I channel!"
        @assert oL < ver4.loopNum "LVer loopNum must be smaller than the ver4 loopNum"

        lvl = ver4.level + 1
        oR = ver4.loopNum - 1 - oL
        TinLidx = ver4.Tidx
        TinRidx = ver4.Tidx + (oL + 1)
        # println("create left/right vertex4 for channel $chan, loopNum $(ver4.loopNum)")

        if chan == T
            if ver4.inBox == false
                Lver = _Ver4(lvl, oL, F, TinLidx, LEFT, false)
                Rver = _Ver4(lvl, oR, FULL, TinRidx, RIGHT, false)
            else
                Lver = _Ver4(lvl, oL, F_CT, TinLidx, LEFT, true)
                Rver = _Ver4(lvl, oR, FULL_CT, TinRidx, RIGHT, true)
            end
        elseif chan == U
            @assert ver4.inBox == false "Ver4 in box can't have U diagrams!"
            Lver = _Ver4(lvl, oL, F, TinLidx, LEFT, false)
            Rver = _Ver4(lvl, oR, FULL, TinRidx, RIGHT, false)
        elseif chan == S
            @assert ver4.inBox == false "Ver4 in box can't have S diagrams!"
            Lver = _Ver4(lvl, oL, V, TinLidx, LEFT, false)
            Rver = _Ver4(lvl, oR, FULL, TinRidx, RIGHT, false)
        elseif chan == TC || chan == UC
            Lver = _Ver4(lvl, oL, F_CT, TinLidx, LEFT, true)
            Rver = _Ver4(lvl, oR, FULL_CT, TinRidx, RIGHT, true)
        end

        # println("para: ", TinLidx, ", ", Lver.Tidx, ", ", ver4.inBox)
        @assert Lver.Tidx == ver4.Tidx "Lver Tidx must be equal to vertex4 Tidx! LoopNum: $(ver4.loopNum), LverLoopNum: $(Lver.loopNum), chan: $chan"

        ############## construct IdxMap ########################################
        map = []
        for (lt, LvT) in enumerate(Lver.Tpair)
            for (rt, RvT) in enumerate(Rver.Tpair)
                GTidx = addTidx(ver4.G[1], (LvT[OUTR], RvT[INL]))
                GTxidx, VerTidx = (0, 0)

                if chan == T
                    VerTidx = addTidx(ver4, (LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR]))
                    GTxidx = addTidx(ver4.G[chan], (RvT[OUTL], LvT[INR]))
                elseif chan == U
                    VerTidx = addTidx(ver4, (LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL]))
                    GTxidx = addTidx(ver4.G[chan], (RvT[OUTL], LvT[INR]))
                elseif chan == S
                    VerTidx = addTidx(ver4, (LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR]))
                    GTxidx = addTidx(ver4.G[chan], (LvT[OUTL], RvT[INR]))
                elseif chan == TC || chan == UC
                    # counterterms are equal-time
                    VerTidx = addTidx(ver4, (LvT[INL], LvT[INL], LvT[INL], LvT[INL]))
                    GTxidx = addTidx(ver4.G[chan], (RvT[OUTL], LvT[INR]))
                end

                @assert ver4.Tpair[end][1] == ver4.Tidx "InL Tidx must be the same for all Tpairs in the vertex4"
                @assert GTidx != 0 && VerTidx != 0 "index cann't be zero!"
                push!(map, IdxMap(lt, rt, GTidx, GTxidx, VerTidx))
            end
        end
        return new(chan, Lver, Rver, map)
    end
end

struct Ver4
    level::Int
    loopNum::Int
    chan::Vector{Int} # list of channels
    Tidx::Int # inital Tidx
    side::Int # right side vertex is always a full gamma4
    inBox::Bool

    G::Vector{Green}
    Tpair::Vector{Tuple{Int,Int,Int,Int}}
    weight::Vector{VerWeight}
    bubble::Vector{Bubble{Ver4}}

    function Ver4(lvl, loopNum, chan, tidx, side, inbox, isfast = false)
        g = [Green() for i = 1:6]
        ver4 = new(lvl, loopNum, chan, tidx, side, inbox, g, [], [], [])
        if loopNum <= 0
            # negative loopNum should never be used in evaluation
            addTidx(ver4, (tidx, tidx, tidx, tidx))
            push!(ver4.weight, zero(VerWeight))
            return ver4
        end
        UST = [c for c in chan if c != I]
        II = [c for c in chan if c == I]
        for c in UST
            for ol = 0:loopNum - 1
                bubble = Bubble{Ver4}(ver4, c, ol)
                if length(bubble.map) > 0
                    push!(ver4.bubble, bubble)
                end
            end
        end
        Num = (isfast && lvl == 0) ? 4 : length(ver4.Tpair)
        for t = 1:Num
            push!(ver4.weight, zero(VerWeight))
        end
        for g in ver4.G
            for t in g.Tpair
                push!(g.weight, 0.0)
            end
        end
        return ver4
    end
end

function eval(
    ver4::Ver4,
    KinL::Mom,
    KoutL::Mom,
    KinR::Mom,
    KoutR::Mom,
    Kidx::Int,
    fast = false,
)
    if ver4.loopNum == 0
        DiagType == POLAR ?
        ver4.weight[1] = interaction(KinL, KoutL, KinR, KoutR, ver4.inBox, norm(varK[0])) :
        ver4.weight[1] = interaction(KinL, KoutL, KinR, KoutR, ver4.inBox)
        return
    end

    # LoopNum>=1
    ver4.weight .*= 0.0 # initialize all weights
    G = ver4.G
    eval(G[1], varK[Kidx])

    for c in ver4.chan
        if c == T || c == TC
            ver4.inBox == false && eval(G[T], KoutL + G[1].K - KinL)
        elseif c == U || c == UC
            # can not be in box!
            eval(G[U], KoutR + G[1].K - KinL)
        else
            # S channel, and cann't be in box!
            eval(G[S], KinL + KinR - G[1].K)
        end
    end

    w = zero(VerWeight)
    gWeight, projfactor = (0.0, 0.0)
    for b in ver4.bubble
        c = b.chan
        projfactor = SymFactor[c] * PhaseFactor

        if c == T || c == TC
            eval(b.Lver, KinL, KoutL, G[T].K, G[1].K, Kidx + 1)
            eval(b.Rver, G[1].K, G[T].K, KinR, KoutR, Kidx + 1)
        elseif c == U || c == UC
            eval(b.Lver, KinL, KoutR, G[U].K, G[1].K, Kidx + 1)
            eval(b.Rver, G[1].K, G[U].K, KinR, KoutL, Kidx + 1)
        else
            # S channel
            eval(b.Lver, KinL, G[S].K, KinR, G[1].K, Kidx + 1)
            eval(b.Rver, G[1].K, KoutL, G[S].K, KoutR, Kidx + 1)
        end


        for map in b.map
            gWeight = (c == TC || c == UC || ver4.inBox) ? counterBubble(G[1].K) :
                G[1].weight[map.G] * G[c].weight[map.Gx]

            Lw, Rw = (b.Lver.weight[map.Lver], b.Rver.weight[map.Rver])

            if c == T || c == TC
                w[DIR] = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR]
                w[EX] = Lw[EX] * Rw[EX]
            elseif c == U || c == UC
                w[DIR] = Lw[EX] * Rw[EX]
                w[EX] = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR]
            else
                # S channel,  see the note "code convention"
                w[DIR] = Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR]
                w[EX] = Lw[DIR] * Rw[DIR] + Lw[EX] * Rw[EX]
            end

            weight = ver4.weight

            if fast && ver4.level == 0
                pair = ver4.Tpair[map.ver]
                dT = varT[pair[INL]] - varT[pair[OUTL]] + varT[pair[INR]] - varT[pair[OUTR]]
                gWeight *= projfactor * cos(2.0 * pi / Beta * dT)
                weight[ChanMap[c]] .+= w .* gWeight
            else
                weight[map.ver] .+= w .* (gWeight * projfactor)
            end
        end
    end
end

function _expandBubble(children, text, style, bub::Bubble, parent)
    push!(children, zeros(Int, 0))
    @assert parent == length(children)
    dict = Dict(
        I => ("I", "yellow"),
        T => ("T", "red"),
        TC => ("Tc", "pink"),
        U => ("U", "blue"),
        UC => ("Uc", "navy"),
        S => ("S", "green"),
    )
    push!(text, "$(dict[bub.chan][1])\n$(bub.Lver.loopNum)-$(bub.Rver.loopNum)")
    push!(style, "fill:$(dict[bub.chan][2])")

    current = length(children) + 1
    push!(children[parent], current)
    _expandVer4(children, text, style, bub.Lver, current) # left vertex 

    current = length(children) + 1
    push!(children[parent], current)
    _expandVer4(children, text, style, bub.Rver, current) # right vertex
end

function _expandVer4(children, text, style, ver4::Ver4, parent)
    push!(children, zeros(Int, 0))
    @assert parent == length(children)
    # println("Ver4: $(ver4.level), Bubble: $(length(ver4.bubble))")
    if ver4.loopNum > 0
        info = "O$(ver4.loopNum)\nT[$(length(ver4.Tpair))]\n"
        # for t in ver4.Tpair
        #     info *= "[$(t[1]) $(t[2]) $(t[3]) $(t[4])]\n"
        # end
    else
        info = "O$(ver4.loopNum)"
    end
    push!(text, info)

    ver4.inBox ? push!(style, "stroke-dasharray:3,2") : push!(style, "")

    for bub in ver4.bubble
        # println(bub.chan)
        current = length(children) + 1
        push!(children[parent], current)
        _expandBubble(children, text, style, bub, current)
    end
end

function visualize(ver4::Ver4)
    children, text, style = (Vector{Vector{Int}}(undef, 0), [], [])
    _expandVer4(children, text, style, ver4, 1)

    # text = ["one\n(second line)", "2", "III", "four"]
    # style = ["", "fill:red", "r:14", "opacity:0.7"]
    # link_style = ["", "stroke:blue", "", "stroke-width:10px"]
    tooltip = ["pops", "up", "on", "hover"]
    t = D3Trees.D3Tree(
        children,
        text = text,
        style = style,
        tooltip = tooltip,
        # link_style = link_style,
        title = "Vertex4 Tree",
        init_expand = 2,
    )

    D3Trees.inchrome(t)
end

export init, eval, Ver4, Bubble, visualize, addTidx

end
