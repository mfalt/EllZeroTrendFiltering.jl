
module dev

import Base.-
import Base.+
import Base.show
import Base.Operators
import Base.start
import Base.next
import Base.done
import Base.length

import Base.==

import IterTools

using StaticArrays
using PyPlot
using Polynomials
using QuadGK


include("types/QuadraticPolynomial.jl")
include("types/PiecewiseQuadratic.jl")
include("types/QuadraticForm2.jl")




function roots(p::QuadraticPolynomial)
    quad_expr = p.b^2 - 4*p.a*p.c
    if quad_expr < 0
        return NaN, NaN
    end

    if p.a != 0
        term1 = (p.b / 2 / p.a)
        term2 = sqrt(quad_expr) / 2 / abs(p.a)
        return -term1-term2, -term1+term2
    else
        term = -p.c / p.b
        return term, Inf
    end
end






# TODO come up with better symbol for ρ
"""
Inserts a quadratic polynomial ρ into the linked list Λ which represents a piecewise quadratic polynomial
"""
function add_quadratic{T}(Λ::PiecewiseQuadratic{T}, ρ::QuadraticPolynomial{T})

    if isnull(Λ.next)
        insert(Λ, ρ, -1e9)
        return
    end


    λ_prev = Λ
    λ_curr = Λ.next

    Δ = QuadraticPolynomial(0., 0., 0.)
    while ~isnull(λ_curr)

        left_endpoint = get(λ_curr).left_endpoint

        # TODO: This is probably not needed now..
        if left_endpoint == -1e9
            #println("minus-inf")
            left_endpoint = -10000.0
        end

        right_endpoint = get_right_endpoint(get(λ_curr))

        if right_endpoint == 1e9
            #println("inf")
            right_endpoint = left_endpoint + 20000.0
        end

        Δ .= ρ .- get(λ_curr).p

        root1,root2 = roots(Δ)

        #println(root1, " : ", root2, " .... (", (left_endpoint, right_endpoint))

        if root1 > left_endpoint && root1 < right_endpoint
            #println("case 1:")
            λ_prev, λ_curr = impose_quadratic_to_root(λ_prev, get(λ_curr), root1, (left_endpoint+root1)/2, ρ, Δ)
            #println(Λ)
        end

        if isnull(λ_curr)
            break
        end

        if root2 > left_endpoint && root2 < right_endpoint
            #println("case 2:")
            λ_prev, λ_curr = impose_quadratic_to_root(λ_prev, get(λ_curr), root2, (root1 + root2)/2, ρ, Δ)
            #println(Λ)
            if Δ.a > 0; return; end # Saves perhaps 5% of computation time
        end

        if isnull(λ_curr)
            break
        end

        #println("case 3:")
        λ_prev, λ_curr = impose_quadratic_to_endpoint(λ_prev, get(λ_curr), (get(λ_curr).left_endpoint + right_endpoint)/2, ρ, Δ)
        #println(Λ)

        if isnull(λ_curr)
            break
        end
    end

end



# TODO Possibly just insert all roots, and then do the comparisons?
# Leads to some unnecessary nodes being created but is perhaps simplier?
function impose_quadratic_to_root{T}(λ_prev::PiecewiseQuadratic{T}, λ_curr::PiecewiseQuadratic{T}, root::T, midpoint::T, ρ::QuadraticPolynomial{T}, Δ::QuadraticPolynomial{T})
    if Δ((λ_curr.left_endpoint + root)/2) < 0 # ρ is smallest, i.e., should be inserted
        if λ_prev.p === ρ
            #println("1")
            λ_curr.left_endpoint = root
            return λ_prev, Nullable{PiecewiseQuadratic{Float64}}(λ_curr)
        else
            #println("2")

            #λ_new = insert(λ_prev, ρ, λ_curr.left_endpoint)
            #λ_curr.left_endpoint = root
            #return λ_new, λ_curr

            λ_new = insert(λ_curr, λ_curr.p, root)
            λ_curr.p = ρ
            return λ_curr, λ_new
        end
    else
        #println("3")
        λ_new = insert(λ_curr, λ_curr.p, root) # insert new interval, but defer deciding on wheather λ_curr.p or ρ is smallest
        return λ_curr, λ_new
    end
end

@inline function impose_quadratic_to_endpoint{T}(λ_prev::PiecewiseQuadratic{T}, λ_curr::PiecewiseQuadratic{T}, midpoint, ρ, Δ)
    local v1, v2
    if Δ(midpoint) < 0  # ρ is smallest, i.e., should be inserted
        if λ_prev.p === ρ
            #println("1")
            λ_new = delete_next(λ_prev)
            v1, v2 = λ_prev, λ_new
        else
            #println("2")
            λ_curr.p = ρ
            v1, v2 = λ_curr, λ_curr.next
        end
    else
        # Do nothing
        #println("3")
        v1, v2 = λ_curr, λ_curr.next
    end
    return v1, v2 #Single point of exit
end

function minimize_wrt_x2(qf::QuadraticForm2)
    P = qf.P
    q = qf.q
    r = qf.r

    if P[2,2] > 0
        QuadraticPolynomial(P[1,1] - P[1,2]^2 / P[2,2],
        q[1] - P[1,2]*q[2] / P[2,2],
        r - q[2]^2 / P[2,2]/ 4)
    elseif P[2,2] == 0 || P[1,2] == 0 || q[2] == 0
        QuadraticPolynomial(P[1,1], q[1], r)
    else
        # There are some special cases, but disregards these
        QuadraticPolynomial(0.,0.,-Inf)
    end
end

function minimize_wrt_x2_fast{T}(qf::QuadraticForm2{T},a,b,c)
    P = qf.P
    q = qf.q
    r = qf.r

    if P[2,2] > 0
        QuadraticPolynomial(P[1,1] - P[1,2]^2 / (P[2,2]+a),
        q[1] - P[1,2]*(q[2]+b) / (P[2,2]+a),
        (r+c) - (q[2]+b)^2 / (P[2,2]+a)/ 4)
    elseif (P[2,2]+a) == 0 || P[1,2] == 0 || (q[2]+b) == 0
        QuadraticPolynomial(P[1,1], q[1], r+c)
    else
        # There are some special cases, but disregards these
        QuadraticPolynomial(0.,0.,-Inf)
    end
end


"""
Find optimal fit
"""
function find_optimal_fit{T}(Λ_0::Array{PiecewiseQuadratic{T},1}, ℓ::Array{QuadraticForm2{T},2}, M)

    N = size(ℓ, 2)

    Λ = Array{PiecewiseQuadratic{T}}(M, N-1)

    Λ[1, :] .= Λ_0

    for m=2:M
        for i=N-m+1:-1:1
            Λ_new = create_new_pwq()
            for ip=i+1:N-m

                for λ in Λ[m-1, ip]
                    p = λ.p

                    # ρ = dev.minimize_wrt_x2(
                    # ℓ[i,ip] + dev.QuadraticForm2{T}(@SMatrix([0. 0; 0 1])*p.a, @SVector([0., 1])*p.b, p.c))
                    # # Avoid ceting two extra QuadraticForm2
                    ρ = dev.minimize_wrt_x2_fast(ℓ[i,ip], p.a, p.b, p.c)
                    ρ.time_index = ip
                    ρ.ancestor = p


                    dev.add_quadratic(Λ_new, ρ)
                end
            end
            Λ[m, i] = Λ_new
        end
    end
    return Λ
end


# Finds the minimum of a positive definite quadratic one variable polynomial
# the find_minimum fcn returns (opt_x, opt_val)
function find_minimum(p::QuadraticPolynomial)
    if p.a <= 0
        println("No unique minimum exists")
        return (NaN, NaN)
    else
        x_opt = -p.b / 2 / p.a
        f_opt = -p.b^2/4/p.a + p.c
        return (x_opt, f_opt)
    end
end


function find_minimum(Λ::PiecewiseQuadratic)
    # May assume that the minimum is at the staionary point of the polynomial
    # TODO: True?

    f_opt = Inf
    p_opt = Λ.p # The segment containing the smallest polynimal
    x_opt = NaN

    for λ in Λ
        x, f = find_minimum(λ.p)
        if f < f_opt
            f_opt = f
            x_opt = x
            p_opt = λ.p
        end
    end
    return p_opt, x_opt, f_opt
end


function recover_solution{T}(Λ::PiecewiseQuadratic{T}, first_index=1, last_index=-1)

    p, y, f = find_minimum(Λ)

    I = [first_index]

    while true
        push!(I, p.time_index)

        if isnull(p.ancestor);
            break;
        end

        p = get(p.ancestor)
    end

    if last_index != -1
        I[end] = last_index
    end

    return I, y, f
end



# TODO Maybe use big T for time indices, howabout mathcal{T}
"""
Computes the transition costs ℓ given a
polynomial and a sequence t
"""
function compute_transition_costs(g, t::AbstractArray)
    T = Float64
    # Find primitive functions to g, t*g, and g^2
    # and evaluate them at the break points
    #I_g = polyint(g).(t)
    #I_g2 = polyint(g^2).(t)
    #I_tg = polyint(Poly([0,1]) * g).(t)

    N = length(t)

    # Find primitive functions to g, t*g, and g^2 at the break points
    I_g = zeros(size(t))
    I_g2 = zeros(size(t))
    I_tg = zeros(size(t))

    for i=2:N
        I_g[i] = I_g[i-1] + quadgk(g, t[i-1], t[i])[1]
        I_g2[i] = I_g2[i-1] + quadgk(t -> g(t)^2, t[i-1], t[i])[1]
        I_tg[i] = I_tg[i-1] + quadgk(t -> t*g(t), t[i-1], t[i])[1]
    end


    ℓ = Array{QuadraticForm2{T}}(N-1,N)

    for i=1:N-1
        for ip=i+1:N

            P = (t[ip] - t[i]) * @SMatrix [1/3 1/6; 1/6 1/3]

            q = -2* 1/(t[ip]-t[i]) *
            @SVector [-(I_tg[ip] - I_tg[i]) + t[ip]*(I_g[ip] - I_g[i]),
            (I_tg[ip] - I_tg[i]) - t[i]*(I_g[ip] - I_g[i])]

            r =  I_g2[ip] - I_g2[i]

            ℓ[i,ip] = QuadraticForm2(P, q, r)
        end
    end

    return ℓ
end

"""
Evaluate the optimal cost (using least squares) for all
possible index sets with K elements
"""
function brute_force_optimization(ℓ, K)
    cost_best = Inf

    I_best = []
    Y_best = []

    N = size(ℓ, 2)

    for I=IterTools.subsets(2:N-1, K)

        I = [1; I; N]
        P = zeros(K+2, K+2)
        q = zeros(K+2)
        r = 0

        # Form quadratic cost function Y'*P*Y + q'*Y + r
        # corresponding to the y-values in the vector Y
        for j=1:K+1
            P[j:j+1,j:j+1] .+= ℓ[I[j], I[j+1]].P
            q[j:j+1] .+= ℓ[I[j], I[j+1]].q
            r += ℓ[I[j], I[j+1]].r
        end

        # find the optimal Y-vector, and compute the correspinding error
        Yopt = -(P \ q) / 2
        cost = Yopt' * P * Yopt + q' * Yopt + r

        if cost < cost_best
            cost_best = cost
            Y_best = Yopt
            I_best = I
        end
    end
    return I_best, Y_best, cost_best
end

function find_optimal_y_values(ℓ, I)
    P = zeros(length(I), length(I))
    q = zeros(length(I))
    r = 0

    # Form quadratic cost function Y'*P*Y + q'*Y + r
    # corresponding to the y-values in the vector Y
    for j=1:length(I)-1
        P[j:j+1,j:j+1] .+= ℓ[I[j], I[j+1]].P
        q[j:j+1] .+= ℓ[I[j], I[j+1]].q
        r += ℓ[I[j], I[j+1]].r
    end

    # find the optimal Y-vector, and compute the correspinding error
    Y = -(P \ q) / 2
    f = Y' * P * Y + q' * Y + r
    return Y, f
end


# end of module
end
