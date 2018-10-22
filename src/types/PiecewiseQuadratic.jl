export create_new_pwq, generate_PiecewiseQuadratic, insert, delete_next
export get_right_endpoint, evalPwq, get_vals, find_minimum, find_minimum_value

"""
    PiecewiseQuadratic{T}

Represents a piecewise-quadratic as a linked list of QuadraticPolynomial{T} and left endpoints.

Note that the first element of the list is sometimes assumed to be a list head with `NaN` polynomial and `left_endpoint=NaN`.

The list/element `x` is the last element in the list when `x.next === x` or `x.next.left_endpoint==Inf`.
"""
mutable struct PiecewiseQuadratic{T}
    π::QuadraticPolynomial{T}
    left_endpoint::T
    next::PiecewiseQuadratic{T}

    PiecewiseQuadratic{T}()  where T = (x=new(); x.left_endpoint=Inf; x.next = x)
    PiecewiseQuadratic{T}(p, left_endpoint, next) where T = new(p, left_endpoint, next)
end

"""
Constructs piece of quadratic polynomial poiting to NULL, i.e.
a rightmost segment, or one that has not been inserted into a piecewise quadratic
"""
function PiecewiseQuadratic(p::QuadraticPolynomial{T}, left_endpoint) where T
    PiecewiseQuadratic{T}(p, left_endpoint, PiecewiseQuadratic{T}())
end
function PiecewiseQuadratic(p::QuadraticPolynomial{T}, left_endpoint, next::PiecewiseQuadratic{T}) where T
    PiecewiseQuadratic{T}(p, left_endpoint, next)
end

islast(λ::PiecewiseQuadratic) = λ.next.left_endpoint == Inf

"""
Creates empty piecewise-quadratic polynomial consisting only of the list head
"""
function create_new_pwq(T)
    return PiecewiseQuadratic(QuadraticPolynomial(T[NaN,NaN,NaN]), T(NaN), PiecewiseQuadratic{T}())
end

"""
Creates piecewise-quadratic polynomial containing one element
"""
function create_new_pwq(p::QuadraticPolynomial{T}) where T
    pwq = PiecewiseQuadratic(p, T(-Inf), PiecewiseQuadratic{T}())
    return PiecewiseQuadratic(QuadraticPolynomial(T[NaN,NaN,NaN]), T(NaN), pwq)
end

function generate_PiecewiseQuadratic(args::Tuple{Vector{T},T}...) where T
    return PiecewiseQuadratic(QuadraticPolynomial(T[NaN,NaN,NaN]), T(NaN), _generate_PiecewiseQuadratic_helper(args...))
end


function _generate_PiecewiseQuadratic_helper(arg)
    return PiecewiseQuadratic(QuadraticPolynomial(arg[1]), arg[2])
end

function _generate_PiecewiseQuadratic_helper(arg, args...)
    return PiecewiseQuadratic(QuadraticPolynomial(arg[1]), arg[2],  _generate_PiecewiseQuadratic_helper(args...))
end

function Base.iterate(pwq::PiecewiseQuadratic{T}, iterstate::PiecewiseQuadratic{T}=pwq.next) where T
    if (iterstate.left_endpoint == Inf)
        return nothing
    end

    return (iterstate, iterstate.next)
end


# For trouble-shooting etc.
function getindex(Λ::PiecewiseQuadratic, n::Integer)
    if n <= 0
        throw(BoundsError(Λ, n))
    end
    # If first element has NaN left_endpoint, assume it to be empty list-head
    λ = isnan(Λ.left_endpoint) ? Λ.next : Λ
    for k=1:n-1
        if islast(λ)
            throw(BoundsError(Λ, n))
        end
        λ = λ.next
    end
    return λ
end


function insert(pwq::PiecewiseQuadratic{T}, p::QuadraticPolynomial, left_endpoint) where T
    pwq.next = PiecewiseQuadratic(p, left_endpoint, pwq.next)
    return pwq.next
end

# Delete the node after pwq and return the node that will follow after
# pwq after the deletion
#OBS This function is unsafe if pwq.next does not exist
function delete_next(pwq::PiecewiseQuadratic{T}) where T
    pwq.next = pwq.next.next
    return pwq.next
end

function get_right_endpoint(λ::PiecewiseQuadratic)
    return λ.next.left_endpoint
end


function length(pwq::PiecewiseQuadratic)
    n = 0
    for x in pwq
        n += 1
    end
    return n
end


function show(io::IO, Λ::PiecewiseQuadratic{T}) where T
    print(io, "PiecewiseQuadratic{$T} with $(length(Λ)) elements:\n")

    for λ in Λ
        if λ.left_endpoint == -Inf
            print(io, "[  -∞ ,")
        else
            @printf(io, "[%3.2f,", λ.left_endpoint)
        end

        if get_right_endpoint(λ) == Inf
            print(io, " ∞  ]")
        else
            @printf(io, " %3.2f]", get_right_endpoint(λ))
        end

        print(io, "\t  :   ")
        show(io, λ.π)
        print(io, "\n")
    end
    return
end

function (Λ::PiecewiseQuadratic)(x::Number)
    for λ in Λ
        if λ.left_endpoint <= x < get_right_endpoint(λ)
            return λ.π(x)
        end
    end
    return NaN
end

function (Λ::PiecewiseQuadratic{T1})(x::AbstractArray) where T1
    y = fill(zero(promote_type(T1,eltype(x))), length(x))
    for λ in Λ
        inds = λ.left_endpoint .<= x .< get_right_endpoint(λ)
        y[inds] .= λ.π.(x[inds])
    end
    return y
end

"""
    x, y, x_all, y_all = get_vals(Λ::PiecewiseQuadratic)
Evaluate `Λ` so that `y .= Λ.(x)` for a set of gridpoints `x`.
Also returns lists x_all, y_all so that `y_all[i] .= Λ[i].(x_all[i])`
where `Λ[i]` is the `i`th segment of the piecewise quadratic.

The gridpoints in x_all[i] cointains points in the interval `[xmin[i]-1,xmax[i]+1]` where
`xmin[i],xmax[i]` are the extrema of the interval for which `Λ[i]` defines `Λ`.

The gridpoints at the ends default to `xmin[1] := xmax[1] - 1` and `xmax[end] := xmin[end] + 1`.
"""
function get_vals(Λ::PiecewiseQuadratic)
    x = Float64[]
    y = Float64[]
    x_all = Array{Float64,1}[]
    y_all = Array{Float64,1}[]
    for λ in Λ
        left_endpoint = λ.left_endpoint
        right_endpoint = get_right_endpoint(λ)

        if left_endpoint != -Inf || right_endpoint != Inf
            if left_endpoint == -Inf
                left_endpoint = right_endpoint - 2.0
            end

            if right_endpoint == Inf
                right_endpoint = left_endpoint + 2.0
            end
        else
            left_endpoint = -10.
            right_endpoint = 10.
        end

        y_grid = range(left_endpoint, stop=right_endpoint, length=50)
        vals =  λ.π.(y_grid)
        append!(x, y_grid)
        append!(y, vals)
        push!(x_all, y_grid)
        push!(y_all, vals)
    end
    return x, y, x_all, y_all
end

#
# function plot_pwq(Λ::PiecewiseQuadratic)
#     x, y, x_all, y_all = get_vals(Λ)
#     p = plot(x,y, l =(3,:red), lab="minimum", size=(1200,800))
#     for i in eachindex(x_all)
#         plot!(p, x_all[i][[1,end]], y_all[i][[1,end]], l = 0, m=(1,:cross,:orange))
#         plot!(p, x_all[i], y_all[i], l=(1,:blue), lab="")
#     end
#     p
#     return p
# end

"""
    p_opt, x_opt, f_opt = find_minimum(Λ::PiecewiseQuadratic)

Find the point `x_opt = argmin Λ(x)` and return `f_opt = Λ(x_opt)` and the quadratic
function `p_opt` which defines `Λ` around the point `x_opt`.
"""
function find_minimum(Λ::PiecewiseQuadratic)
    # May assume that the minimum is at the staionary point of the polynomial
    # TODO: True?

    f_opt = Inf
    π_opt = typeof(Λ.π)()
    x_opt = NaN

    for λ in Λ
        x, f = find_minimum(λ.π)
        if f < f_opt
            f_opt = f
            x_opt = x
            π_opt = λ.π
        end
    end
    return π_opt, x_opt, f_opt
end

function find_minimum_value(Λ::PiecewiseQuadratic)
    f_opt = Inf

    for λ in Λ
        _, f = find_minimum(λ.π)
        if f < f_opt
            f_opt = f
        end
    end
    if isnan(f_opt)
        println("find_minimum_value failed to find a minimum, no minimum exists")
    end
    return f_opt
end
