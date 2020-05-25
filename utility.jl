module Utility
# using StaticArrays
# mutable struct VerWeight <: StaticArray{Tuple{2},Float64,1}
#     dir::Float64
#     ex::Float64
#     VerWeight(x, y) = new(x, y)
# end

# @inline function Main.getindex(v::VerWeight, i::Int)
#     T = eltype(v)
#     return GC.@preserve v unsafe_load(Base.unsafe_convert(Ptr{T}, pointer_from_objref(v)), i)
# end

# @inline function Main.setindex!(v::VerWeight, val, i::Int)
#     T = eltype(v)
#     GC.@preserve v unsafe_store!(Base.unsafe_convert(Ptr{T}, pointer_from_objref(v)), convert(T, val), i)
#     return val
# end

# export VerWeight
end
