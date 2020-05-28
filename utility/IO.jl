import JLD2, FileIO, Statistics, Measurements

function loadFile(folder, filename)
    Norm = []
    Data = []
    for file in readdir(folder)
        # println(file)
        if occursin(filename, file)
            filepath =
            JLD2.jldopen(joinpath(folder, file), "r") do f
                println("Loading ", file)
                push!(Norm, f["Norm"])
                push!(Data, f["Data"])
            end
        end
    end
    return Data, Norm
end

function statis(data, weight)
    N = length(data)
    @assert N == length(weight) "data and weight should have the same length"
    @assert N > 0 "data is empty!"
    Z = sum(weight)
    avg = sum(data .* weight) / Z
    # println(avg)
    var = sum((d - avg).^2 .* w / Z for (d, w) in zip(data, weight))
    # println(var)
    err = sqrt.(var / N)
    return Measurements.measurement.(avg, err)
    # return avg, var
end