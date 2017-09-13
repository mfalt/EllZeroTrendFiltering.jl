
module dev

import Base.-
import Base.+
import Base.show
import Base.Operators
import Base: start, next, done, length, zero, getindex, ==

import IterTools
using StaticArrays
using QuadGK


include("types/QuadraticPolynomial.jl")
include("types/PiecewiseQuadratic.jl")
include("types/QuadraticForm2.jl")

include("transition_cost_computation.jl")

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
function find_optimal_fit{T}(ℓ::Array{QuadraticForm2{T},2}, V_0N::QuadraticPolynomial{T}, M::Integer, upper_bound=Inf)
    #global counter1
    #global counter2
    #counter1 = 0
    #counter2 = 0

    N = size(ℓ, 2)

    Λ = Array{PiecewiseQuadratic{T}}(M, N)

    for i=1:N-1
        p = dev.minimize_wrt_x2(ℓ[i, N], V_0N)
        p.time_index = N
        Λ[1, i] .= dev.create_new_pwq(p)
    end

    ρ = QuadraticPolynomial{T}()
    upper_bound_inner = Inf
    for m=2:M
        for i=1:N-m
            Λ_new = create_new_pwq()
            for ip=i+1:N-m+1

                for λ in Λ[m-1, ip]
                    p = λ.p

                    #counter1 += 1

                    minimize_wrt_x2(ℓ[i,ip], p, ρ)
                    ρmin = unsafe_minimum(ρ)
                    if ρmin > upper_bound || ρmin > upper_bound_inner
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
            if i == 1
                _, _, upper_bound_inner = find_minimum(Λ[m,1])
            end
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


function recover_optimal_index_set{T}(Λ::PiecewiseQuadratic{T}, first_index=1, last_index=-1)

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

"""
Y, f = find_optimal_y_values(ℓ, V_0N, I)
Given transition costs `ℓ`, cost of right endpoint `V_0N`, and breakpoint indicies `I`
the optimal y-values `Y` and the optimal cost `f` are computed.
"""

function find_optimal_y_values(ℓ, V_0N::QuadraticPolynomial, I)
    m = length(I) - 1

    P = zeros(m+1, m+1)
    q = zeros(m+1)

    # Add cost for the right endpoint
    P[m+1, m+1] = V_0N.a
    q[m+1] = V_0N.b
    r = V_0N.c

    # Form quadratic cost function Y'*P*Y + q'*Y + r
    # corresponding to the y-values in the vector Y
    for j=1:m
        P[j:j+1,j:j+1] .+= ℓ[I[j], I[j+1]].P
        q[j:j+1] .+= ℓ[I[j], I[j+1]].q
        r += ℓ[I[j], I[j+1]].r
    end

    # find the optimal Y-vector, and compute the correspinding error
    Y = -(P \ q) / 2
    f = Y' * P * Y + q' * Y + r
    return Y, f
end

"""
    I, Y, f = recover_solution(Λ::PiecewiseQuadratic{T}, ℓ, V_0N::QuadraticPolynomial, first_index=1)

"""
function recover_solution(Λ::PiecewiseQuadratic, ℓ, V_0N::QuadraticPolynomial, first_index=1)
    I, _, _ = recover_optimal_index_set(Λ, first_index)
    Y, f = find_optimal_y_values(ℓ, V_0N::QuadraticPolynomial, I)

    return I, Y, f
end

"""
Evaluate the optimal cost (using least squares) for all
possible index sets with M segemets
"""
function brute_force_optimization(ℓ, V_0N::QuadraticPolynomial, M)
    cost_best = Inf

    I_best = []
    Y_best = []

    N = size(ℓ, 2)

    for I=IterTools.subsets(2:N-1, M-1)

        I = [1; I; N]
        P = zeros(M+1, M+1)
        q = zeros(M+1)
        r = 0

        # Add cost at right endpoint 
        P[end,end] = V_0N.a
        q[end]     = V_0N.b
        r          = V_0N.c
        # Form quadratic cost function Y'*P*Y + q'*Y + r
        # corresponding to the y-values in the vector Y
        for j=1:M
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




# end of module
end
