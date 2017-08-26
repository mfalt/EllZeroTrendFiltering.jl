# Quadratic Polynomial to represent cost-to-go function,
# keeos track of its "ancestor" to simplify the recovery of the solution
# TODO: penality for making mutable, so ancestor can be set later on,
# or set ancestor at creation
mutable struct QuadraticPolynomial{T<:Real}
a::T
b::T
c::T
ancestor::Nullable{QuadraticPolynomial{T}}
time_index::Int
function QuadraticPolynomial{T}(a::T, b::T, c::T) where {T}
    # @assert a ≥ 0, # Δ may have negative a ...
    new(a, b, c, Nullable{QuadraticPolynomial{T}}(),-1)
end
function QuadraticPolynomial{T}(a::T, b::T, c::T,ancestor, time_index) where {T}
    @assert a ≥ 0
    new(a, b, c, ancestor, time_index)
end
end

QuadraticPolynomial{T}(a::T,b::T,c::T) = QuadraticPolynomial{T}(a,b,c)

function QuadraticPolynomial(x::Vector)
    QuadraticPolynomial(x[3], x[2], x[1])
end

function Base.show(io::IO, p::QuadraticPolynomial)
    print(@sprintf("%.2f*x^2 + %.2f*x + %.2f   ", p.a, p.b, p.c))
end




-(p1::QuadraticPolynomial, p2::QuadraticPolynomial) = QuadraticPolynomial(p1.a-p2.a, p1.b-p2.b, p1.c-p2.c)


==(p1::QuadraticPolynomial, p2::QuadraticPolynomial) = (p1.a==p2.a) && (p1.b==p2.b) && (p1.c==p2.c)


polyval(p::QuadraticPolynomial, x) = p.a*x^2 + p.b*x + p.c

function (p::QuadraticPolynomial)(x)
    return polyval(p, x)
end



# Finds the minimum of a positive definite quadratic one variable polynomial
# the find_minimum fcn returns (opt_x, opt_val)
function find_minimum(p::QuadraticPolynomial)
    if p.a <= 0
        println("No unique minimum exists")
        return (NaN, NaN)
    else
        return (-p.b / 2 / p.a, -p.b^2/4 + p.c)
    end
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
