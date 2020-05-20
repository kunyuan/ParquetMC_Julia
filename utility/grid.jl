# asserting() = false # when set to true, this will enable all `@myassert`s

# macro tmpassert(test)
#     esc(:(if $(@__MODULE__).asserting()
#         @assert($test)
#     end))
# end

module Grid
using StaticArrays

struct Coeff 
    bound::SVector{2,Real}
    idx::SVector{2,Real}
    lambda::Real
    a::Real
    b::Real

    function Coeff(b, i, l, dense2sparse)
        bound = b
        idx = i
        lambda = l

        if dense2sparse == false
            bound[1], bound[2] = bound[2], bound[1]
            idx[1], idx[2] = idx[2], idx[1]
            lambda *= -1
        end
        _l1, _l2 = 1.0, exp(lambda * (idx[2] - idx[1]))
        b = (bound[2] - bound[1]) / (_l2 - _l1)
        a = (bound[1] * _l2 - bound[2] * _l1) / (_l2 - _l1)
        return new(bound, idx, lambda, a, b)
    end
end

function _floor(l::Coeff, x::Real)::Int
    # @tmpassert(bound[0] <= x <= bound[1])
    pos = l.idx[1] + 1.0 / l.lambda * log((x - l.a) / l.b)
    return Base.floor(Int, pos)
end

function _grid(l::Coeff, idx::Int)::Real
    return l.a + l.b * exp(l.lambda * (Real(idx) - l.idx[1]))
end

struct LogGrid
    size::Int
    grid::Array{Real,1}
    segment::Array{Int,1}
    coeff::Array{Coeff,1}

    function LogGrid(coeff::Array{Coeff,1}, segment::Array{Int,1})
        @assert segment[1] == 1 "Segment shoudl start with 1!"
        @assert length(segment) == length(coeff) + 1 "Number of coeff and segments don't match!"
        grid = Real[]
        for i in 1:length(segment) - 1
            for idx in segment[i]:segment[i + 1]
                push!(grid, _grid(coeff[i], idx))
            end
        end
        return new(length(grid), grid, segment, coeff)
    end
end

struct UniformGrid
    size::Int
    grid::Array{Real,1}
    
    function UniformGrid(head::Real, tail::Real, size::Int)
        @assert size > 1 "Size must be large than 1"
        grid = LinRange(head, tail, size)
        return new(size, grid)
    end
end

    export Coeff, LogGrid, UniformGrid
end

# using .Grid

# logGrid = Grid.LogGrid([0.0 10.0], [1.0, 32.0], 2.0 / 32.5, true)
# for i in 1:32
#     println(Grid.grid(logGrid, i))
# end

