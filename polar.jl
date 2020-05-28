include("parameter.jl")
include("utility/IO.jl")
include("grid.jl")
using PyCall
# import Trapz:trapz
pushfirst!(PyVector(pyimport("sys")."path"), "")
pushfirst!(PyVector(pyimport("sys")."path"), joinpath(pwd(), "utility/"))

Data, Norm = loadFile("Data", r"polar_pid[0-9]+")
F = pyimport("fourier")
plt = pyimport("matplotlib.pyplot")

Fourier = F.fourier(Grid.tau.grid, [0.0, ], Beta)
MomGrid = Grid.K.grid
TauGrid = Grid.tau.grid
println(size(Data[1][1, 1, :]))

for o in 1:Order
    yList = [Fourier.naiveT2W(PyCall.PyReverseDims(d[:, :, o])) for d in Data]
    # yList = [Fourier.naiveT2W(d[:, :, o]) for d in Data]
    polar, err = statis(yList, Norm)
    # polar, err = [], []
    # for (idx, k) in enumerate(MomGrid)
    #     data = [trapz(TauGrid, d[1, idx, :]) for d in Data]
    #     y, e = statis(data, Norm)
    #     push!(polar, y)
    #     push!(err, e)
    # end
    plt.errorbar(MomGrid / Kf, polar, yerr = err, fmt = "o-", capthick = 1, capsize = 4,
                    label = "Order $o")
end

# y, err = statis(Data, Norm)
# println(polar[1, :, 1])
# using PyPlot
# plt = PyPlot.plot(Grid.K.grid, y[1, :, 1])
# xlabel!("q")
# display(plt)
# gui()
# readline()

# plt.plot(Grid.K.grid, y[1, :, 1])
plt.legend()
plt.show()


# println(Norm)

