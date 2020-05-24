module Vertex4
include("../parameter.jl")
include("propagator.jl")
using .Propagator: interaction, green, counterBubble
using D3Trees

const ChanMap = [I, T, U, S, T, U]
const SymFactor = [1.0, -1.0, 1.0, -0.5, 1.0, -1.0]

struct Green
    Tpair::Vector{Tuple{Int,Int}}
    weight::Vector{Float}
    Green() = new([], [])
end

function eval(G::Green, K::Mom, _varT, IsAnomal = false)
    for (i, t) in enumerate(G.Tpair)
        G.weight[i] = green(_varT[t[OUT]] - _varT[t[IN]], K)
    end
end

function addTidx(obj, _Tidx)
    for (i, Tidx) in enumerate(obj.Tpair)
        if Tidx == _Tidx
            return i
        end
    end
    push!(obj.Tpair, _Tidx)
    return length(obj.Tpair)
end

struct IdxMap
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
        G = ver4.G
        for (lt, LvT) in enumerate(Lver.Tpair)
            for (rt, RvT) in enumerate(Rver.Tpair)
                GT0idx = addTidx(G[1], (LvT[OUTR], RvT[INL]))
                GTxidx, VerTidx = (0, 0)

                if chan == T
                    VerTidx = addTidx(ver4, (LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR]))
                    GTxidx = addTidx(G[T], (RvT[OUTL], LvT[INR]))
                elseif chan == U
                    VerTidx = addTidx(ver4, (LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL]))
                    GTxidx = addTidx(G[U], (RvT[OUTL], LvT[INR]))
                elseif chan == S
                    VerTidx = addTidx(ver4, (LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR]))
                    GTxidx = addTidx(G[S], (LvT[OUTL], RvT[INR]))
                elseif chan == TC || chan == UC
                    # counterterms are equal-time
                    VerTidx = addTidx(ver4, (LvT[INL], LvT[INL], LvT[INL], LvT[INL]))
                    GTxidx = -1 # do not need it
                # GTxidx = addTidx(ver4.G[chan], (RvT[OUTL], LvT[INR]))
                else
                    throw("This channel is invalid!")
                end

                @assert ver4.Tpair[end][1] == ver4.Tidx "InL Tidx must be the same for all Tpairs in the vertex4"
                @assert GT0idx != 0 && VerTidx != 0 "index cann't be zero!"
                push!(map, IdxMap(GT0idx, GTxidx, VerTidx))
            end
        end
        return new(chan, Lver, Rver, map)
    end
end

struct Ver4
    level::Int
    loopNum::Int
    chan::Set{Int} # list of channels
    Tidx::Int # inital Tidx
    side::Int # right side vertex is always a full gamma4
    inBox::Bool

    K::MVector{3,Mom}
    G::MVector{4,Green}
    Tpair::Vector{Tuple{Int,Int,Int,Int}}
    weight::Vector{VerWeight}
    bubble::Vector{Bubble{Ver4}}

    function Ver4(lvl, loopNum, chan, tidx, side, inbox, isfast = false)
        g = @MVector [Green() for i = 1:4]
        k = @MVector [zero(Mom) for i = 1:3]
        ver4 = new(lvl, loopNum, Set(chan), tidx, side, inbox, k, g, [], [], [])
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

function eval(ver4::Ver4, KinL, KoutL, KinR, KoutR, Kidx::Int, var, fast = false)
    if ver4.loopNum == 0
        DiagType == POLAR ?
        ver4.weight[1] = interaction(KinL, KoutL, KinR, KoutR, ver4.inBox, norm(var.K[0])) :
        ver4.weight[1] = interaction(KinL, KoutL, KinR, KoutR, ver4.inBox)
        return
    end

    # LoopNum>=1
    ver4.weight .*= 0.0 # initialize all weights
    G = ver4.G
    K, Kt, Ku, Ks = (var.K[Kidx], ver4.K[1], ver4.K[2], ver4.K[3])
    eval(G[1], K, var.T)
    bubWeight = counterBubble(K)

    for c in ver4.chan
        if c == T || c == TC
            Kt .= KoutL .+ K .- KinL
            if (!ver4.inBox)
                eval(G[T], Kt, var.T)
            end
        elseif c == U || c == UC
            # can not be in box!
            Ku .= KoutR .+ K .- KinL
            eval(G[U], Ku, var.T)
        else
            # S channel, and cann't be in box!
            Ks .= KinL .+ KinR .- K
            eval(G[S], Ks, var.T)
        end
    end

    for b in ver4.bubble
        c = b.chan
        Factor = SymFactor[c] * PhaseFactor

        if c == T || c == TC
            eval(b.Lver, KinL, KoutL, Kt, K, Kidx + 1, var)
            eval(b.Rver, K, Kt, KinR, KoutR, Kidx + 1, var)
        elseif c == U || c == UC
            eval(b.Lver, KinL, KoutR, Ku, K, Kidx + 1, var)
            eval(b.Rver, K, Ku, KinR, KoutL, Kidx + 1, var)
        else
            # S channel
            eval(b.Lver, KinL, Ks, KinR, K, Kidx + 1, var)
            eval(b.Rver, K, KoutL, Ks, KoutR, Kidx + 1, var)
        end

        rN = length(b.Rver.weight)
        for (l, Lw) in enumerate(b.Lver.weight)
            for (r, Rw) in enumerate(b.Rver.weight)
                map = b.map[(l - 1) * rN + r]

                if ver4.inBox || c == TC || c == UC
                    gWeight = bubWeight * Factor
                else
                    gWeight = G[1].weight[map.G] * G[c].weight[map.Gx] * Factor
                end

                if fast && ver4.level == 0
                    pair = ver4.Tpair[map.ver]
                    dT =
                        var.T[pair[INL]] - var.T[pair[OUTL]] + var.T[pair[INR]] -
                        var.T[pair[OUTR]]
                    gWeight *= cos(2.0 * pi / Beta * dT)
                    w = ver4.weight[ChanMap[c]]
                else
                    w = ver4.weight[map.ver]
                end

                if c == T || c == TC
                    w.dir +=
                        gWeight * (Lw.dir * Rw.dir * SPIN + Lw.dir * Rw.ex + Lw.ex * Rw.dir)
                    w.ex += gWeight * Lw.ex * Rw.ex
                elseif c == U || c == UC
                    w.dir += gWeight * Lw.ex * Rw.ex
                    w.ex +=
                        gWeight * (Lw.dir * Rw.dir * SPIN + Lw.dir * Rw.ex + Lw.ex * Rw.dir)
                else
                    # S channel,  see the note "code convention"
                    w.dir += gWeight * (Lw.dir * Rw.ex + Lw.ex * Rw.dir)
                    w.ex += gWeight * (Lw.dir * Rw.dir + Lw.ex * Rw.ex)
                end

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
        for t in ver4.Tpair
            info *= "[$(t[1]) $(t[2]) $(t[3]) $(t[4])]\n"
        end
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

export Ver4, Bubble, visualize

end
