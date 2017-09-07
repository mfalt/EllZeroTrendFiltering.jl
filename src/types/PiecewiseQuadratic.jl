# Note that the first element of the list should be interpreted as the list head
# Should perhaps be reconsidered
#
# Also note that left and right endpoints are represented with ±1e9,
# in order to allow comparison of polynomials far to the left/right,
# Could also just look at the sign of Δ.a
mutable struct PiecewiseQuadratic{T}
p::QuadraticPolynomial{T}
left_endpoint::Float64
next::PiecewiseQuadratic{T}
PiecewiseQuadratic{T}()  where T = (x=new(); x.left_endpoint=Inf; x.next = x)
PiecewiseQuadratic{T}(p, left_endpoint, next) where T = new(p, left_endpoint, next)
end

"""
Constructs piece of quadratic polynomial poiting to NULL, i.e.
a rightmost piece, or one that has not been inserted into a piecewise quadfatic
"""
function PiecewiseQuadratic{T}(p::QuadraticPolynomial{T}, left_endpoint)
    PiecewiseQuadratic{T}(p, left_endpoint, PiecewiseQuadratic{T}())
end
function PiecewiseQuadratic{T}(p::QuadraticPolynomial{T}, left_endpoint, next::PiecewiseQuadratic{T})
    PiecewiseQuadratic{T}(p, left_endpoint, next)
end


"""
Creates empty piecewise-quadratic polynomial consisting only of the list head
"""
function create_new_pwq()
    return PiecewiseQuadratic(QuadraticPolynomial([NaN,NaN,NaN]), NaN, PiecewiseQuadratic{Float64}())
end

"""
Creates piecewise-quadratic polynomial containing one element
"""
function create_new_pwq(p::QuadraticPolynomial)
    pwq = PiecewiseQuadratic(p, -1e9, PiecewiseQuadratic{Float64}())
    return PiecewiseQuadratic(QuadraticPolynomial([NaN,NaN,NaN]), NaN, pwq)
end

function generate_PiecewiseQuadratic(args)
    return PiecewiseQuadratic(QuadraticPolynomial([NaN,NaN,NaN]), NaN, _generate_PiecewiseQuadratic_helper(args))
end

function _generate_PiecewiseQuadratic_helper(args)
    if Base.length(args) > 1
        return PiecewiseQuadratic(QuadraticPolynomial(args[1][1]), args[1][2],  _generate_PiecewiseQuadratic_helper(args[2:end]))
    else
        return PiecewiseQuadratic(QuadraticPolynomial(args[1][1]), args[1][2])
    end
end


start{T}(pwq::PiecewiseQuadratic{T}) = pwq.next
done{T}(pwq::PiecewiseQuadratic{T}, iterstate::PiecewiseQuadratic{T}) = (iterstate.left_endpoint == Inf)
next{T}(pwq::PiecewiseQuadratic{T}, iterstate::PiecewiseQuadratic{T}) = (iterstate, iterstate.next)


# For trouble-shooting etc.
function getindex(Λ::dev.PiecewiseQuadratic, n::Int64)
    if n <= 0
        error("Attempted to access index <= 0 in piecewise quadratic.")
    end

    λ = get(Λ.next)
    for k=1:n-1
        λ = get(λ.next)
    end
    return λ
end


function insert{T}(pwq::PiecewiseQuadratic{T}, p::QuadraticPolynomial, left_endpoint)
    pwq.next = PiecewiseQuadratic(p, left_endpoint, pwq.next)
    return pwq.next
end

# Delete the node after pwq and return the node that will follow after
# pwq after the deletion
#OBS This function is unsafe if pwq.next does not exist
function delete_next{T}(pwq::PiecewiseQuadratic{T})
    pwq.next = unsafe_get(pwq.next).next
    return pwq.next
end

function get_right_endpoint(λ::PiecewiseQuadratic)
    if isnull(λ.next)
        return 1e9
    else
        return unsafe_get(λ.next).left_endpoint
    end
end


function length(pwq::PiecewiseQuadratic)
    n = 0
    for x in pwq # Skip the dummy head element
        n += 1
    end
    return n
end


function show(io::IO, Λ::PiecewiseQuadratic)
    #@printf("PiecewiseQuadratic (%i elems):\n", length(Λ))

    for λ in Λ
        if λ.left_endpoint == -1e9
            print(io, "[  -∞ ,")
        else
            @printf(io, "[%3.2f,", λ.left_endpoint)
        end

        if get_right_endpoint(λ) == 1e9
            print(io, " ∞  ]")
        else
            @printf(io, " %3.2f]", get_right_endpoint(λ))
        end

        print(io, "\t  :   ")
        show(io, λ.p)
        print(io, "\n")
    end
    return
end


function evalPwq(Λ::PiecewiseQuadratic, x)
    y = zeros(x)
    for λ in Λ
        inds = λ.left_endpoint .<= x .< get_right_endpoint(λ)
        y[inds] .= λ.p.(x[inds])
    end
    return y
end

function get_vals(Λ::PiecewiseQuadratic)
    x = Float64[]
    y = Float64[]
    x_all = Array{Float64,1}[]
    y_all = Array{Float64,1}[]
    for λ in Λ
        left_endpoint = λ.left_endpoint
        right_endpoint = get_right_endpoint(λ)

        if left_endpoint != -1e9 || right_endpoint != 1e9
            if left_endpoint == -1e9
                left_endpoint = right_endpoint - 2.0
            end

            if right_endpoint == 1e9
                right_endpoint = left_endpoint + 2.0
            end
        else
            left_endpoint = -10.
            right_endpoint = 10.
        end

        y_grid_gray = linspace(left_endpoint-1, right_endpoint+1)
        push!(x_all, y_grid_gray)
        push!(y_all, λ.p.(y_grid_gray))

        y_grid = linspace(left_endpoint, right_endpoint)
        append!(x, y_grid)
        append!(y, λ.p.(y_grid))
    end
    return x, y, x_all, y_all
end


function plot_pwq(Λ::PiecewiseQuadratic)
    x, y, x_all, y_all = get_vals(Λ)
    p = plot(x,y, l =(3,:red), lab="minimum", size=(1200,800))
    for i in eachindex(x_all)
        plot!(p, x_all[i][[1,end]], y_all[i][[1,end]], l = 0, m=(1,:cross,:orange))
        plot!(p, x_all[i], y_all[i], l=(1,:blue), lab="")
    end
    p
    return p
end






"""
min!(q1::PiecewiseQuadratic{T}, q2::QuadraticOnInterval{T})
Update `q1(x)` on the interval `I` in `q2` to be `q1(x) := min(q1(x),q2(x)) ∀ x∈I`
"""
