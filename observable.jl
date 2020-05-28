
module Observable
include("parameter.jl")
include("grid.jl")
using StaticArrays: MArray
using JLD2, FileIO

const curr = Main.Curr

mutable struct OneBody
    norm::Float # normalization 
    phy::Float # the physcial weight of the normalization diagrams
    # static::MArray{(Order, KGridSize),Float}
    data::Array{Float,3} # order, kgrid, taugrid
    OneBody() = new(1.0e-10, 1.0, zeros(Float, Order, KGridSize, TauGridSize))
end

function measure(obs::OneBody, weight, factor)
    if curr.order == 0
        obs.norm += weight * factor
    else
        obs.data[curr.order, curr.extKidx, curr.extTidx] += weight * factor
    end
end

function save(obs::OneBody)
    filename = "$(name())_pid$(curr.PID).jld2"
    data = Dict("PID" => curr.PID, "Norm" => obs.norm, "Data" => obs.data / obs.norm * obs.phy)

    # FileIO.save(filename, data, compress = true)
    FileIO.save(filename, data)
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
