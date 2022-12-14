export PowerIO

abstract type LinearIO end

struct PowerIO <: LinearIO
    xVec::Vector{Float64}
    yVec::Vector{Float64}
    wVec::Vector{Float64}
    XMat::Matrix{Float64}
    AMat::Matrix{Float64}

    function PowerIO(xdat::VecI{Tx}, ydat::VecI{Ty}, n::Int) where {Tx<:Real, Ty<:Real}
        xVec = Vector{Float64}(undef, n)
        yVec = Vector{Float64}(undef, n)
        wVec = Vector{Float64}(undef, 2)
        XMat = Matrix{Float64}(undef, 2, 2)
        AMat = Matrix{Float64}(undef, 2, 2)

        copyto!(xVec, xdat)
        copyto!(yVec, ydat)
        @inbounds begin
            AMat[1,1] = n
            AMat[1,2] = AMat[2,1] = sum(xdat)
            AMat[2,2] = dot(xdat, xdat, n)

            detA = AMat[1,1] * AMat[2,2] - AMat[1,2] * AMat[2,1]
            AMat[1,1], AMat[2,2] = AMat[2,2] / detA, AMat[1,1] / detA
            AMat[1,2] /= -detA
            AMat[2,1] /= -detA

            wVec1 = sum(ydat)
            wVec2 = dot(xdat, ydat, n)
            wVec[1] = AMat[1,1] * wVec1 + AMat[1,2] * wVec2
            wVec[2] = AMat[2,1] * wVec1 + AMat[2,2] * wVec2
        end

        new(xVec, yVec, wVec, XMat, AMat)
    end
end

function PowerIO(; xdat::VecI=nothing, ydat::VecI=nothing)
    n = length(xdat)
    n ≠ length(ydat) && error("PowerIO: length(xdat) ≠ length(ydat)")
    return PowerIO(xdat, ydat, n)
end

@inline (o::PowerIO)(x::Real) = _transform(o.wVec, x)
@inline _transform(w::VecI, x::Real) = @inbounds w[1] + w[2] * x
