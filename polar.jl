include("parameter.jl")
include("utility/IO.jl")
using Measurements

Data, Norm = loadFile("Data", r"polar_pid[0-9]+")
println(size(Data[1]))
Data = [sum(d, dims = 3) for d in Data]
polar = statis(Data, Norm)
println(polar[1, :, 1])
# println(Norm)

