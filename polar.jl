include("parameter.jl")
include("utility/IO.jl")
include("grid.jl")
# using PyCall
import Trapz:trapz
using Plots
# pushfirst!(PyVector(pyimport("sys")."path"), "")
# pushfirst!(PyVector(pyimport("sys")."path"), joinpath(pwd(), "utility/"))

Data, Norm = loadFile("Data", r"polar_pid[0-9]+")
# F = pyimport("fourier")
# plt = pyimport("matplotlib.pyplot")

# Fourier = F.fourier(Grid.tau.grid, [0.0, ], Beta)
MomGrid = Grid.K.grid
TauGrid = Grid.tau.grid
# println(MomGrid)
# println(size(Data[1][1, 1, :]))
println(Data[1][:, 1, 1])
# println(Data[1][2, :, 1])
println(length(TauGrid))

# for ki in 1:KGridSize
#     println("k=$(Grid.K.grid[ki]): ", sum(Data[1][:, ki, 1]) / length(TauGridSize) * Beta / obs.norm * obs.phy)
# end

# for o in 1:Order
    # yList = [Fourier.naiveT2W(PyCall.PyReverseDims(d[:, :, o])) for d in Data]
    # yList = [Fourier.naiveT2W(d[:, :, o]) for d in Data]
    # polar, err = statis(yList, Norm)
    # size(polar)
    polar, err = [], []
    
    for (idx, k) in enumerate(MomGrid)
        data = [trapz(TauGrid, d[:, idx, 1]) for d in Data]
        # data = [sum(d[:, idx, 1]) / length(TauGrid) * Beta for d in Data]
        # data = [sum(d[:, idx, 1])  * Beta for d in Data]
        # println(data[1])
        y, e = statis(data, Norm)
        push!(polar, y)
        push!(err, e)
    end

    # plot(MomGrid/Kf, polar)
    # plt.errorbar(MomGrid / Kf, polar, yerr = 0.0, fmt = "o-", capthick = 1, capsize = 4, label = "Order $o")
    # y = sum(Data[1][:, :, 1], dims = 1) * Beta / length(TauGrid)
    # println(length(y))
    # println(size(y))
    # plt.errorbar(MomGrid / Kf, y[1, :], yerr = 0.0, fmt = "o-", capthick = 1, capsize = 4, label = "Order $o")

    # y, e = statis([d[:, 20, 1] for d in Data], Norm)
    # plt.errorbar(TauGrid, y, yerr = e, fmt = "o-", capthick = 1, capsize = 4,
    #                 label = "Order $o")
    # y, e = statis([d[1, :, 1] for d in Data], Norm)
    # plt.errorbar(MomGrid / Kf, y, yerr = e, fmt = "o-", capthick = 1, capsize = 4, label = "Order $o")
# end
plot(1:10, 1:10)
display()


# for k in 1:4:length(MomGrid)
#     data = [d[:, k, 1] for d in Data]
#     y, err = statis(data, Norm)
#     plt.errorbar(TauGrid, y, yerr = err, fmt = "o-", capthick = 1, capsize = 4,
#                     label = "K $(MomGrid[k] / Kf)")
# end

# y, err = statis(Data, Norm)
# println(polar[1, :, 1])
# using PyPlot
# plt = PyPlot.plot(Grid.K.grid, y[1, :, 1])
# xlabel!("q")
# display(plt)
# gui()
# readline()

# plt.plot(Grid.K.grid, y[1, :, 1])
# plt.legend()
# plt.show()


# println(Norm)

