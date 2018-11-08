# Quadratic form of 2 variables, used for representing the transition costs l
struct QuadraticForm{T}
    P::SMatrix{2,2,T,4}
    q::SVector{2,T}
    r::T
end
#TODO Remove when StaticArrays properly handle scalars
function QuadraticForm(P::SMatrix{2,2,T,4}, q::AbstractVecOrMat{T}, r::T) where T
    return QuadraticForm(P, SVector{2,T}(q), r)
end

function find_minimum(qf::QuadraticForm)
    x_opt = -qf.P \ pf.q / 2
    f_opt = -qf.q⋅(qf.P\qf.q)/4+qf.r
    return x_opt, f_opt
end

(qf::QuadraticForm)(y::Number, yp::Number) = [y,yp]'*qf.P*[y,yp] + qf.q'*[y, yp] + qf.r
(qf::QuadraticForm)(x::AbstractVector) = x'*qf.P*x + qf.q'*x + qf.r


+(qf1::QuadraticForm, qf2::QuadraticForm) = QuadraticForm(qf1.P+qf2.P, qf1.q+qf2.q, qf1.r+qf2.r)
==(qf1::QuadraticForm, qf2::QuadraticForm) = (qf1.P==qf2.P) && (qf1.q==qf2.q) && (qf1.r==qf2.r)


"""
    μ = minimize_wrt_x2(qf::QuadraticForm{T}, p::QuadraticPolynomial{T}, μ=QuadraticPolynomial{T}()) where {T}
Takes a quadratic form in `[x₁; x₂]` and a polynomial in `x₂`
and returns the minimum of the sum wrt to `x₂`,
i.e. μ(x₁) = min_x₂{ qf(x₁,x₂) +  p(x₂) }`

The input `μ` can be pre-allocated on input and will then be changed.
"""
@inline function minimize_wrt_x2_old(qf::QuadraticForm{T}, p::QuadraticPolynomial{T}, μ=QuadraticPolynomial()) where T
    P = qf.P
    q = qf.q
    r = qf.r

    P22_new = P[2,2] + p.a

    if P22_new > 0
        μ.a = P[1,1] - P[1,2]^2 / P22_new
        μ.b = q[1] - P[1,2]*(q[2] + p.b) / P22_new
        μ.c = (r + p.c) - (q[2] + p.b)^2 / P22_new / 4
    elseif P22_new == 0
        μ.a = P[1,1]
        μ.b = q[1]
        μ.c = r + p.c
    else
        error("Piecewise quadratic cost should be positive (semi-)definite. Please submit a bug report.")
    end
    return μ
end


@inline function minimize_wrt_x2_alt(qf::QuadraticForm{T}, χ::SVector{2,T}, p::QuadraticPolynomial{T}, μ=QuadraticPolynomial{T}()) where T
    P = qf.P + p.a*χ*χ'
    q = qf.q + p.b*χ
    r = qf.r + p.c

    if P[2,2] > 0
        μ.a = P[1,1] - P[1,2]^2 / P[2,2]
        μ.b = q[1] - P[1,2]*q[2] / P[2,2]
        μ.c = r - q[2]^2 / P[2,2] / 4
    elseif P[2,2] == 0
        μ.a = P[1,1]
        μ.b = q[1]
        μ.c = r
    else
        error("Piecewise quadratic cost should be positive (semi-)definite. Please submit a bug report.")
    end
    return μ
end


@inline function minimize_wrt_x2(qf::QuadraticForm{T}, μ=QuadraticPolynomial{T}()) where T
    P = qf.P
    q = qf.q
    r = qf.r

    if P[2,2] > 0
        μ.a = P[1,1] - P[1,2]^2 / P[2,2]
        μ.b = q[1] - P[1,2]*q[2] / P[2,2]
        μ.c = r - q[2]^2 / P[2,2] / 4
    elseif P[2,2] == 0
        μ.a = P[1,1]
        μ.b = q[1]
        μ.c = r
    else
        error("Piecewise quadratic cost should be positive (semi-)definite. Please submit a bug report.")
    end
    return μ
end
