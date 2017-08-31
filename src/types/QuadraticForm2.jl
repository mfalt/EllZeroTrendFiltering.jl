# Quadratic form of 2 variables, used for representing the transition costs ℓ
type QuadraticForm2{T}
    P::SMatrix{2,2,T,4}
    q::SVector{2,T}
    r::T
end

function (qf::QuadraticForm2)(y, yp)
    return [y,yp]'*qf.P*[y,yp] + qf.q'*[y, yp] + qf.r
end


+(qf1::QuadraticForm2, qf2::QuadraticForm2) = QuadraticForm2(qf1.P+qf2.P, qf1.q+qf2.q, qf1.r+qf2.r)
==(qf1::QuadraticForm2, qf2::QuadraticForm2) = (qf1.P==qf2.P) && (qf1.q==qf2.q) && (qf1.r==qf2.r)
