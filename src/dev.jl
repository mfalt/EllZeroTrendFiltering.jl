
module dev

import Base.-
import Base.+
import Base.show
import Base.Operators
import Base.start
import Base.next
import Base.done
import Base.length
import Base.getindex

import Base.==

import IterTools
using StaticArrays
using QuadGK


include("types/QuadraticPolynomial.jl")
include("types/PiecewiseQuadratic.jl")
include("types/QuadraticForm2.jl")

global const DEBUG = false
global const DEBUG2 = false
global const COUNTER_TEST = false

# TODO come up with better symbol for ρ
"""
Inserts a quadratic polynomial ρ into the linked list Λ which represents a piecewise quadratic polynomial
"""
function add_quadratic{T}(Λ::PiecewiseQuadratic{T}, ρ::QuadraticPolynomial{T})


    DEBUG && println("Inserting: ", ρ)
    if Λ.next.left_endpoint == Inf # I.e. the piecewise quadratic object is empty, perhaps better to add dummy polynomial
        ρ.has_been_used = true
        insert(Λ, ρ, -Inf)
        return
    end

    λ_prev = Λ
    λ_curr = Λ.next

    while λ_curr.left_endpoint != Inf #???
        #global counter2 += 1
        DEBUG && println(Λ)

        left_endpoint = λ_curr.left_endpoint
        right_endpoint = get_right_endpoint(λ_curr)


        Δa = ρ.a - λ_curr.p.a
        Δb = ρ.b - λ_curr.p.b
        Δc = ρ.c - λ_curr.p.c

        b2_minus_4ac =  Δb^2 - 4*Δa*Δc

        if Δa > 0 # ρ has greater curvature, i.e., ρ is smallest in the middle
            if b2_minus_4ac <= 0
                # No intersections, old quadratic is smallest, just step forward
                DEBUG && println("No intersections, old quadratic is smallest, Δa > 0, breaking.")
                break
            else

                # Compute the intersections
                term1 = -(Δb / 2 / Δa)
                term2 = sqrt(b2_minus_4ac) / 2 / Δa # Δa > 0
                root1, root2 = term1-term2, term1+term2

                DEBUG && println("Δa > 0   root1:", root1, "   root2:", root2)

                # Check where the intersections are and act accordingly
                if root1 >= right_endpoint
                    DEBUG && println("Two intersections to the right")
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                elseif root2 <= left_endpoint
                    # No intersections, old quadratic is smallest, step forward
                    DEBUG && println("Two intersections to the left")
                    break # There will be no more intersections since Δa > 0
                elseif root1 <= left_endpoint && root2 >= right_endpoint
                    # No intersections, new quadratic is smallest
                    DEBUG && println("One intersections on either side")
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, ρ)
                elseif root1 > left_endpoint && root2 < right_endpoint
                    DEBUG && println("Two intersections within the interval")
                    λ_prev, λ_curr = update_segment_old_new_old(λ_curr, ρ, root1, root2)
                    break # There will be no more intersections since Δa > 0
                elseif root1 > left_endpoint
                    DEBUG && println("Root 1 within the interval")
                    λ_prev, λ_curr = update_segment_old_new(λ_curr, ρ, root1)
                elseif root2 < right_endpoint
                    DEBUG && println("Root 2 within the interval")
                    λ_prev, λ_curr = update_segment_new_old(λ_prev, λ_curr, ρ, root2)
                    break # There will be no more intersections since Δa > 0
                else
                    error("Shouldn't end up here")
                end
            end

        elseif Δa < 0 # ρ has lower curvature, i.e., ρ is smallest on the sides
            if b2_minus_4ac <= 0
                λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, ρ)
            else
                # Compute the intersections
                term1 = -(Δb / 2 / Δa)
                term2 = sqrt(b2_minus_4ac) / 2 / Δa # < 0
                root1, root2 = term1+term2, term1-term2
                DEBUG && println("Δa < 0   root1:", root1, "   root2:", root2)

                # Check where the intersections are and act accordingly
                if root1 >= right_endpoint || root2 <= left_endpoint
                    # No intersections, ρ is smallest
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, ρ)
                elseif root1 <= left_endpoint && root2 >= right_endpoint
                    # No intersections, old quadratic is smallest, just step forward
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                elseif root1 > left_endpoint && root2 < right_endpoint
                    # Two intersections within the interval
                    λ_prev, λ_curr = update_segment_new_old_new(λ_prev, λ_curr, ρ, root1, root2)
                elseif root1 > left_endpoint
                    λ_prev, λ_curr = update_segment_new_old(λ_prev, λ_curr, ρ, root1)
                elseif root2 < right_endpoint
                    λ_prev, λ_curr = update_segment_old_new(λ_curr, ρ, root2)
                else
                    error("Shouldn't end up here")
                end
            end
        else # a == 0.0
            DEBUG && pritnln("Δa == 0")
            DEBUG2 && println("Δa == 0 : $ρ")
            if Δb == 0
                if Δc >= 0
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                else
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, ρ)
                end
                continue
            end

            root = -Δc / Δb
            if Δb > 0
                if root < left_endpoint
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                elseif root > right_endpoint
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, ρ)
                else
                    λ_prev, λ_curr = update_segment_new_old(λ_prev, λ_curr, ρ, root)
                end
            else
                if root < left_endpoint
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, ρ)
                elseif root > right_endpoint
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                else
                    DEBUG2 && println("Special case")
                    λ_prev, λ_curr = update_segment_old_new(λ_curr, ρ, root)
                end
            end
        end

    end
    return
end



@inline function update_segment_do_nothing(λ_curr)
    return λ_curr, λ_curr.next
end

@inline function update_segment_new(λ_prev, λ_curr, ρ)
    ρ.has_been_used = true
    if λ_prev.p === ρ
        λ_prev.next = λ_curr.next
        v1, v2 = λ_prev, λ_curr.next
    else
        λ_curr.p = ρ
        v1, v2 = λ_curr, λ_curr.next
    end
    return v1, v2
end

@inline function update_segment_old_new(λ_curr, ρ, break1)
    ρ.has_been_used = true
    new_pwq_segment =  PiecewiseQuadratic(ρ, break1, λ_curr.next)
    λ_curr.next = new_pwq_segment
    return new_pwq_segment, new_pwq_segment.next
end

@inline function update_segment_new_old(λ_prev, λ_curr, ρ, break1)
    if λ_prev.p === ρ
        λ_curr.left_endpoint = break1
    else
        ρ.has_been_used = true
        λ_prev.next = PiecewiseQuadratic(ρ, λ_curr.left_endpoint, λ_curr)
        λ_curr.left_endpoint = break1
    end
    return λ_curr, λ_curr.next
end

@inline function update_segment_new_old_new(λ_prev, λ_curr, ρ, break1, break2)
    ρ.has_been_used = true
    update_segment_new_old(λ_prev, λ_curr, ρ, break1)
    return update_segment_old_new(λ_curr, ρ, break2)
end

@inline function update_segment_old_new_old(λ_curr, ρ, break1, break2)
    ρ.has_been_used = true
    second_old_pwq_segment =  PiecewiseQuadratic(λ_curr.p, break2, λ_curr.next)
    new_pwq_segment =  PiecewiseQuadratic(ρ, break1, second_old_pwq_segment)
    λ_curr.next = new_pwq_segment
    return second_old_pwq_segment, second_old_pwq_segment.next
end

# Takes a quadratic form in [x1; x2] and a polynomial in x2
# and returns the minimum of the sum wrt to x2,
# i.e. a polynomial of x1
@inline function minimize_wrt_x2{T}(qf::QuadraticForm2{T},p::QuadraticPolynomial{T},ρ=QuadraticPolynomial{T}())
    P = qf.P
    q = qf.q
    r = qf.r

    P22_new = P[2,2] + p.a

    if P22_new > 0
        ρ.a = P[1,1] - P[1,2]^2 / P22_new
        ρ.b = q[1] - P[1,2]*(q[2]+p.b) / P22_new
        ρ.c = (r+p.c) - (q[2]+p.b)^2 / P22_new / 4
    elseif P22_new == 0 #|| qf.P11 == 0 || (qf.q2+p.b) == 0 #why are the two last conditions needed?
        ρ.a = P[1,1]
        ρ.b = q[1]
        ρ.c = r+p.c
    else
        # FIXME: what are these condtions?
        # There are some special cases, but disregards these
        ρ.a = 0.0
        ρ.b = 0.0
        ρ.c = 0.0
    end
    return ρ
end




global counter1
global counter2

"""
Find optimal fit
"""
function find_optimal_fit{T}(Λ_0::Array{PiecewiseQuadratic{T},1}, ℓ::Array{QuadraticForm2{T},2}, M::Int, upper_bound=Inf)
    #global counter1
    #global counter2
    #counter1 = 0
    #counter2 = 0

    N = size(ℓ, 2)

    Λ = Array{PiecewiseQuadratic{T}}(M, N)

    Λ[1, 1:end-1] .= Λ_0


    ρ = QuadraticPolynomial{T}()
    for m=2:M
        for i=N-m:-1:1
            Λ_new = create_new_pwq()
            for ip=i+1:N-m+1

                for λ in Λ[m-1, ip]
                    p = λ.p

                    #counter1 += 1

                    minimize_wrt_x2(ℓ[i,ip], p, ρ)

                    if unsafe_minimum(ρ) > upper_bound
                        continue
                    end

                    add_quadratic(Λ_new, ρ)

                    if ρ.has_been_used == true
                        ρ.time_index = ip
                        ρ.ancestor = p
                        ρ = QuadraticPolynomial{T}()
                        ρ.has_been_used = false
                    end
                end
            end
            Λ[m, i] = Λ_new
        end
    end
    return Λ
end





"""
Find optimal fit
"""
function regularize{T}(Λ_0::PiecewiseQuadratic{T}, ℓ::Array{QuadraticForm2{T},2}, reg_param::T, upper_bound=Inf)
    N = size(ℓ, 2)

    Λ = Vector{PiecewiseQuadratic{T}}(N)

    Λ[N] = Λ_0

    ρ = QuadraticPolynomial{T}()

    for i=N-1:-1:1
        Λ_new = create_new_pwq()
        for ip=i+1:N

            for λ in Λ[ip]
                p = λ.p
                #counter1 += 1

                minimize_wrt_x2(ℓ[i,ip], p, ρ)
                ρ.c += reg_param

                if unsafe_minimum(ρ) > upper_bound
                    continue
                end

                add_quadratic(Λ_new, ρ)

                if ρ.has_been_used == true
                    ρ.time_index = ip
                    ρ.ancestor = p
                    ρ = QuadraticPolynomial{T}()
                    ρ.has_been_used = false
                end
            end
        end
        Λ[i] = Λ_new
    end

    return Λ
end


function recover_solution{T}(Λ::PiecewiseQuadratic{T}, first_index=1, last_index=-1)

    p, y, f = find_minimum(Λ)

    I = [first_index]

    while true
        push!(I, p.time_index)

        if !isdefined(p, :ancestor);
            break;
        end

        p = p.ancestor
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
        const tol = 1e-3
        I_g[i] = I_g[i-1] + quadgk(g, t[i-1], t[i], reltol=tol)[1]
        I_g2[i] = I_g2[i-1] + quadgk(t -> g(t)^2, t[i-1], t[i], reltol=tol)[1]
        I_tg[i] = I_tg[i-1] + quadgk(t -> t*g(t), t[i-1], t[i], reltol=tol)[1]
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



# Med samma P-matriser?
function compute_discrete_transition_costs(g)
    T = Float64

    N = length(g)

    # Find sums of g, k*g, and g^2
    G1 = zeros(T, N)
    G2 = zeros(T, N)
    G3 = zeros(T, N)

    # The sums corresponding to transitioning from i to ip
    # i.e. not including the cost at ip
    for k=2:N
        G1[k] = G1[k-1] + g[k-1]
        G2[k] = G2[k-1] + (k-1)*g[k-1]
        G3[k] = G3[k-1] + g[k-1]^2
    end

    # The P-matrices only depend on the distance d=ip-i
    P_mats  = Vector{SMatrix{2,2,Float64,4}}(N-1)
    P_mats[1] = @SMatrix [1.0 0; 0 0]
    for d=2:N-1
        off_diag_elems = sum([k*(d - k) for k=0:d-1])
        P_mats[d] = @SMatrix [P_mats[d-1][1,1] + d^2    off_diag_elems;
        off_diag_elems            P_mats[d-1][1,1]]
    end

    P_mats = P_mats ./ (1.0:N-1).^2 # FIXME: Why can't this be done above in the loop?

    #P_invs = inv.(P_mats)

    ℓ = Array{QuadraticForm2{T}}(N-1,N)

    for i=1:N-1
        for ip=i+1:N

            P = P_mats[ip-i]

            q = -2* 1/(ip-i) *
            @SVector [-(G2[ip] - G2[i]) + ip*(G1[ip] - G1[i]),
            (G2[ip] - G2[i]) - i*(G1[ip] - G1[i])]

            r =  G3[ip] - G3[i]

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
