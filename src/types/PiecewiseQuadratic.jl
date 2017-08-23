struct QuadraticOnInterval{T}
    p::QuadraticPolynomial{T}
    endpoint::T
end

struct PiecewiseQuadratic{T}
    qlist::LList{QuadraticOnInterval{T}}
end

PiecewiseQuadratic{T}(p::QuadraticPolynomial{T}) = PiecewiseQuadratic{T}(llist(QuadraticOnInterval(p,typemax(T))))

"""
    min!(q1::PiecewiseQuadratic{T}, q2::QuadraticOnInterval{T})
Update `q1(x)` on the interval `I` in `q2` to be `q1(x) := min(q1(x),q2(x)) ∀ x∈I`
"""
function min!{T}(q1::PiecewiseQuadratic{T}, q2::QuadraticOnInterval{T})
    #TODO
end


Λs = Array{PiecewiseQuadratic,1}(len1,len2)
N = 10
M = 4
for m = M:-1:1
    for i = (N-m):-1:1
        Λim = PiecewiseQuadratic(QuadraticPolynomial(Inf,Inf,Inf))
        for ip = (i+1):(N-m+1)
            update!(Λim, ts[ip], Λ[m,ip], l[i,ip])
