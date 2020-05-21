module Markov
include("parameter.jl")



const UpdateNum = 9

function printStatus()
    ReWeight[2] = 1.6
    println(Beta)
end

function changeOrder()
    return
end

function changeTau()
    return
end

function changeK()
    return
end

function changeExtTau()
    return
end

function changeExtK()
    return
end

function evaluate()
    return
end

export printStatus, changeOrder, changeTau, changeK, changeExtTau, changeExtK, evaluate

end
