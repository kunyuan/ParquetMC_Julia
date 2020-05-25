# include("parameter.jl")
module Test

# import Main.Weight : VerWeight

function init()
    w = Main.VerWeight(0.0, 0.0)
    w.dir = 2.0
    println(typeof(w))
    println(w)
    println(Main.Beta)
end

export init

end