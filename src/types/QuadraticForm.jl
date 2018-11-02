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
