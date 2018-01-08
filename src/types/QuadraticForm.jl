# Quadratic form of 2 variables, used for representing the transition costs ℓ
struct QuadraticForm{T}
    P::SMatrix{2,2,T,4}
    q::SVector{2,T}
    r::T
end

function (qf::QuadraticForm)(y, yp)
    return [y,yp]'*qf.P*[y,yp] + qf.q'*[y, yp] + qf.r
end


+(qf1::QuadraticForm, qf2::QuadraticForm) = QuadraticForm(qf1.P+qf2.P, qf1.q+qf2.q, qf1.r+qf2.r)
==(qf1::QuadraticForm, qf2::QuadraticForm) = (qf1.P==qf2.P) && (qf1.q==qf2.q) && (qf1.r==qf2.r)

function find_minimum_value{T}(qf::QuadraticForm{T})
    if qf.P[2] == qf.P[3] == qf.P[4] == zero(T)
        return -qf.q[1]^2/(4*qf.P[1])+qf.r
    elseif qf.P[1] == qf.P[2] == qf.P[3] == zero(T)
        return -qf.q[2]^2/(4*qf.P[4])+qf.r
    else #Assume P is not zero
        return -qf.q'*(qf.P\qf.q)/4 + qf.r
    end
end
