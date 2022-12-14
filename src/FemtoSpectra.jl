module FemtoSpectra

const VecI = AbstractVector

function dot(x::Vector{Tx}, y::Vector{Ty}, n::Int) where {Tx<:Real, Ty<:Real}
    r = 0.0
    m = mod(n, 5)
    if m ≠ 0
        for i in 1:m
            @inbounds r += x[i] * y[i]
        end
        n < 5 && return r
    end
    m += 1
    for i in m:5:n
        @inbounds r += x[i] * y[i] + x[i+1] * y[i+1] + x[i+2] * y[i+2] + x[i+3] * y[i+3] + x[i+4] * y[i+4]
    end
    return r
end

include("./PowerIO.jl")

end # module
