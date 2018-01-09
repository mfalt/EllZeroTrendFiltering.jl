export find_minimum, roots

# Quadratic Polynomial to represent cost-to-go function,
# keeos track of its "ancestor" to simplify the recovery of the solution
# has_been_used keeps track of if the Quadratic polynimial has been inserted
# into a piecewise quadratic, otherwise it is reused
mutable struct QuadraticPolynomial{T<:Real}
a::T
b::T
c::T
has_been_used::Bool
time_index::Int32
ancestor::QuadraticPolynomial{T}
function QuadraticPolynomial{T}(a::T, b::T, c::T) where {T}
    new(a, b, c, false, -1)
end

QuadraticPolynomial{T}() where {T} = new()
end

QuadraticPolynomial{T}(a::T,b::T,c::T) = QuadraticPolynomial{T}(a,b,c)

function QuadraticPolynomial(x::Vector)
    QuadraticPolynomial(x[3], x[2], x[1])
end

function Base.show(io::IO, p::QuadraticPolynomial)
    @printf(io, "%.2f*x^2 + %.2f*x + %.2f   ", p.a, p.b, p.c)
end

function Base.zero{T}(::Type{QuadraticPolynomial{T}})
    return QuadraticPolynomial(zero(T), zero(T), zero(T))
end
Base.zero{T<:QuadraticPolynomial}(S::T) = zero(T)


-(p1::QuadraticPolynomial, p2::QuadraticPolynomial) = QuadraticPolynomial(p1.a-p2.a, p1.b-p2.b, p1.c-p2.c)

# # x .= y .- z
# # x .= y .+ z
# function Base.broadcast!{T}(op::Function, x::QuadraticPolynomial{T}, y::QuadraticPolynomial{T}, z::QuadraticPolynomial{T})
#     x.a = op(y.a,z.a)
#     x.b = op(y.b,z.b)
#     x.c = op(y.c,z.c)
# end

==(p1::QuadraticPolynomial, p2::QuadraticPolynomial) = (p1.a==p2.a) && (p1.b==p2.b) && (p1.c==p2.c)

@inline function (p::QuadraticPolynomial)(x)
    return p.a*x^2 + p.b*x + p.c
end



# Finds the minimum of a positive definite quadratic one variable polynomial
# the find_minimum fcn returns (opt_x, opt_val)
function find_minimum(p::QuadraticPolynomial)
    if p.a < 0 || (p.a == 0 && p.b != 0)
        println("No unique minimum exists")
        return (NaN, NaN)
    elseif p.a == 0
        # p.b == 0
        return (NaN, p.c)
    else
        x_opt = -p.b / 2 / p.a
        f_opt = -p.b^2/4/p.a + p.c
        return (x_opt, f_opt)
    end
end

"""
Minimum of a quadratic function which is assumed to be positive definite
No checks of this is done
"""
@inline function unsafe_minimum{T}(p::QuadraticPolynomial{T})
    return (-p.b^2/4/p.a + p.c)::T
end
