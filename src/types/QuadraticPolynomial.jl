struct QuadraticPolynomial{T<:Real}
    a::T
    b::T
    c::T
end

function QuadraticPolynomial{T}(a::T,b::T,c::T)
    @assert a ≥ 0
    QuadraticPolynomial{T}(a,b,c)
end

"""
    n, intersec = intersections{T}(p1::QuadraticPolynomial{T},p2::QuadraticPolynomial{T})
    returns number of intersections `n` and a 2-vector `intersec` containing the intersection points
    in the first `n` elements.
"""
function intersections{T}(p1::QuadraticPolynomial{T},p2::QuadraticPolynomial{T})
    a = p1.a-p2.a
    b = p1.b-p2.b
    c = p1.c-p2.c
    n = 0
    intersec = [T(NaN), T(NaN)]
    if a == 0
        #Numerical problems if division, 0 or 1 intersections
        if b == 0
            # 0 intersections
            #return 0, [T(NaN), T(NaN)]
        else
            # 1 intersection
            n = 1
            intersec[1] = -c/b
            #return 1, [-c/b, T(NaN)]
        end
    else
        # 0 or 2 intersections, possibly in same point
        r = b^2 - 4c*a
        if r < 0
            # 0 intersections
            #return 0, [T(NaN), T(NaN)]
        else
            # 2 intersections
            r = sqrt(r)/(2a)                 #sqrt(b^2-4ca)/2a
            v = -b/(2a)
            n = 2
            intersec[1] = v-r
            intersec[2] = v+r
            #return 2, [v-r, v+r]
        end
    end
    return n, intersec
end
