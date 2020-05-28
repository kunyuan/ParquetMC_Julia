
module Observable
include("parameter.jl")
include("grid.jl")
using StaticArrays: MArray

const curr = Main.Curr

mutable struct OneBody
    norm::Float # normalization 
    phy::Float # the physcial weight of the normalization diagrams
    # static::MArray{(Order, KGridSize),Float}
    # estimator::MArray{Tuple{Order,KGridSize,TauGridSize},Float}
    estimator::Array{Float,3} # order, kgrid, taugrid
    function OneBody()
        estimator = zeros(Float, Order, KGridSize, TauGridSize)
        return new(1.0e-10, 1.0, estimator)
    end
end

function measure(obs::OneBody, weight, factor)
    if curr.order == 0
        obs.norm += weight * factor
    else
        obs[curr.order, curr.kidx, curr.tidx] += weight * factor
    end
end

function save(obs::OneBody)
    filename = "$(name())_pid$(curr.PID).dat"
    open(filename, "w") do io
        write(io, "# Counter: $(curr.step)\n")
        write(io, " # Norm: $(obs.norm)\n")
        write(io, " # KGrid: $(Grid.K.grid)\n")
        write(io, " # TauGrid: $Grid.tau.grid)\n")
        for order = 1:Order
            for qidx = 1:KGridSize
                for tidx = 1:TauGridSize
                    write(io, obs.estimator[order, qidx, tidx], " ")
                end
            end
        end
    end
end

function name()
    if DiagType == POLAR
        return "polar"
    elseif DiagType == SIGMA
        return "sigma"
    elseif DiagType == DELTA
        return "delta"
    elseif DiagType == GAMMA
        return "gamma"
    else
        throw("Not implemented!")
    end
end

end
