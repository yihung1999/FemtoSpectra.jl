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

        #=
            ┌              ┐      ┌       ┐
            │ 1⃗ₙ⋅1⃗ₙ  1⃗ₙ⋅x⃗ₙ │      │ 1⃗ₙ⋅y⃗ₙ │
        𝐀 = │              │, 𝐰 = │       │
            │ x⃗ₙ⋅1⃗ₙ  x⃗ₙ⋅x⃗ₙ │      │ x⃗ₙ⋅y⃗ₙ │
            └              ┘      └       ┘
        =#
        @inbounds begin
            XMat[1,1] = n
            XMat[1,2] = XMat[2,1] = sum(xdat)
            XMat[2,2] = dot(xdat, xdat, n)

            detA = XMat[1,1] * XMat[2,2] - XMat[1,2] * XMat[2,1]
            AMat[1,1], AMat[2,2] = XMat[2,2] / detA, XMat[1,1] / detA
            AMat[1,2] = -(XMat[1,2] / detA)
            AMat[2,1] = -(XMat[2,1] / detA)

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

function Base.show(io::IO, o::PowerIO)
    print(io, "PowerIO: A linear regression object\n")
    print(io, "-----------------------------------\n")
    print(io, "  intercept = $(round(o.wVec[1]; sigdigits=6))\n")
    print(io, "  slope = $(round(o.wVec[2]; sigdigits=6))\n")
    return nothing
end

@inline function (o::PowerIO)(x::Real, m::Symbol)
    if m ≡ :F || m ≡ :f || m ≡ :forward
        return _forward_transform(o.wVec, x)
    elseif m ≡ :B || m ≡ :b || m ≡ :backward
        return _backward_transform(o.wVec, x)
    else
        error("PowerIO: $m is an invalid mode.")
    end
end

@inline function (o::PowerIO)(x::Real, m::String)
    m = lowercase(m)
    if m ≡ "f" || m ≡ "forward"
        return _forward_transform(o.wVec, x)
    elseif m ≡ "b" || m ≡ "backward"
        return _backward_transform(o.wVec, x)
    else
        error("PowerIO: $m is an invalid mode.")
    end
end

@inline _forward_transform(w::VecI, x::Real)  = @inbounds w[1] + w[2] * x
@inline _backward_transform(w::VecI, y::Real) = @inbounds (y - w[1]) / w[2]
