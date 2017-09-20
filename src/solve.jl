
# TODO come up with better symbol for ρ
"""
    add_quadratic!(Λ::PiecewiseQuadratic{T}, ρ::QuadraticPolynomial{T}) where {T}
Inserts a quadratic polynomial ρ into the linked list `Λ`, which represents a piecewise
quadratic polynomial, so that the new `Λ` satisfies
`Λ(y) := min{ Λ(y) ,ρ(y) } ∀ y`
"""
function add_quadratic!{T}(Λ::PiecewiseQuadratic{T}, ρ::QuadraticPolynomial{T})


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

        if Δa > 0 # ρ has greater curvature, i.e., ρ is smallest in the middle if intersect
            if b2_minus_4ac <= 0
                # Zero (or one) intersections, old quadratic is smallest, just step forward
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
                # Zero (or one) roots
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
        else # Δa == 0.0
            DEBUG && println("Δa == 0")
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

"""
    ρ = minimize_wrt_x2(qf::QuadraticForm{T}, p::QuadraticPolynomial{T}, ρ=QuadraticPolynomial{T}()) where {T}
Takes a quadratic form in `[x₁; x₂]` and a polynomial in `x₂`
and returns the minimum of the sum wrt to `x₂`,
i.e. ρ(x₁) = min_x₂{ qf(x₁,x₂) +  p(x₂) }`

The input `ρ` can be pre-allocated on input and will then be changed.
"""
@inline function minimize_wrt_x2{T}(qf::QuadraticForm{T}, p::QuadraticPolynomial{T}, ρ=QuadraticPolynomial{T}())
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
    Λ = find_optimal_fit(ℓ::Array{QuadraticForm{T},2}, V_0N::QuadraticPolynomial{T}, M::Integer, upper_bound=Inf) where {T}

Given the transition costs `ℓ[i,j](y_i,y_j)` and the cost at the endpoint `V_0N(y_N)` find all solutions `f` with up to `M` segments for the problem

V_i^m = minimize_f^M [ Σ_{k=1}^i { ℓ[k,k+1](f(k),f(k+1)) } + V_0N(f(N)) ]
s.t          f(k) being continuous piecewise linear with `m` segements.

i.e. `Λ[m,i]` contains the best (in `ℓ` cost) continuous piecewise linear function `f` with up to `M` segments over the interval `i` to `N`
"""
function find_optimal_fit{T}(ℓ::Array{QuadraticForm{T},2}, V_0N::QuadraticPolynomial{T}, M::Integer, upper_bound=Inf)
    #global counter1
    #global counter2
    #counter1 = 0
    #counter2 = 0

    N = size(ℓ, 2)

    @assert M-1 <= N "Cannot have more segments than N-1."

    Λ = Array{PiecewiseQuadratic{T}}(M, N)

    for i=1:N-1
        p = minimize_wrt_x2(ℓ[i, N], V_0N)
        p.time_index = N
        Λ[1, i] .= create_new_pwq(p)
    end

    ρ = QuadraticPolynomial{T}()
    upper_bound_inner = Inf
    global times
    for m=2:M
        #println("m: $m")
        for i=1:N-m
            Λ_new = create_new_pwq()
            if min(upper_bound, upper_bound_inner) < Inf
                OPTIMIZE && add_quadratic!(Λ_new, QuadraticPolynomial{T}(0.0, 0.0, min(upper_bound,upper_bound_inner)))
            end
            #println("m: $m, i: $i")
            for ip=i+1:N-m+1
                DEBUG && println("(m:$m, i:$i, ip:$ip)")
                for λ in Λ[m-1, ip]
                    p = λ.p

                    #counter1 +=

                    minimize_wrt_x2(ℓ[i,ip], p, ρ)

                    DEBUG && println("Obtained ρ = $ρ")

                    ρmin = unsafe_minimum(ρ)
                    if ρmin > upper_bound || ρmin > upper_bound_inner
                        DEBUG && println("Breaking due to that $ρmin > max($upper_bound, $upper_bound_inner)")
                        continue
                    end

                    DEBUG && println("Inserting...")

                    add_quadratic!(Λ_new, ρ)

                    if ρ.has_been_used == true
                        ρ.time_index = ip
                        ρ.ancestor = p
                        ρ = QuadraticPolynomial{T}()
                        ρ.has_been_used = false
                    end
                end
            end
            #remove_over(Λ_new, min(upper_bound_inner,upper_bound))
            Λ[m, i] = Λ_new
            if i == 1
                if OPTIMIZE
                    upper_bound_inner = find_minimum_value(Λ[m,1])
			upper_bound_inner += sqrt(eps())
                end
            end
        end
    end
    return Λ
end

# TODO: regularized
function  fit_pwl_constrained(g::AbstractArray, M)
    ℓ = compute_discrete_transition_costs(g);
    cost_last = QuadraticPolynomial(1.0, -2*g[end], g[end]^2)
    fit_pwl_constrained_internal(ℓ, cost_last, M)
end

function  fit_pwl_constrained(g, t, M)
    ℓ = compute_transition_costs(g, t);
    cost_last = QuadraticPolynomial(0.0, 0.0, 0.0)
    fit_pwl_constrained_internal(ℓ, cost_last, M)
end

function  fit_pwl_constrained_internal(ℓ, cost_last, M)
    Λ = compute_value_functions_constrained(ℓ, cost_last, 0.1);

    Ivec = Vector{Vector{Int}}(M)
    Yvec = Vector{Vector{Float64}}(M)
    fvec = Vector{Float64}(M)

    for m=1:M
        Ivec[m], Yvec[m], fvec[m] = recover_solution(Λ[m, 1], ℓ, cost_last)
    end

    return Ivec, Yvec, fvec
end

"""
Solves the regularization problem
minimzie ∫ (g - y)^2 dt + ζ⋅card(d^2/dt^2 y)
"""
function regularize{T}(ℓ::Array{QuadraticForm{T},2}, V_0N::QuadraticPolynomial{T}, ζ::T)
    N = size(ℓ, 2)

    Λ = Vector{PiecewiseQuadratic{T}}(N)

    V_0N = deepcopy(V_0N)
    V_0N.time_index = -1
    Λ[N] = create_new_pwq(V_0N)

    ρ = QuadraticPolynomial{T}()

    for i=N-1:-1:1
        Λ_new = create_new_pwq()
        for ip=i+1:N


            ζ_level_insertion = false
            for λ in Λ[ip]
                p = λ.p
                #counter1 += 1

                minimize_wrt_x2(ℓ[i,ip], p, ρ)
                ρ.c += ζ # add cost for break point


                add_quadratic!(Λ_new, ρ)


                if ζ_level_insertion == false
                    if !poly_minus_constant_is_greater(Λ_new, ρ, ζ)
                        ζ_level_insertion = true
                    end
                end


                if ρ.has_been_used
                    ρ.time_index = ip
                    ρ.ancestor = p
                    ρ = QuadraticPolynomial{T}()
                    ρ.has_been_used = false

                    ζ_level_insertion = true
                else
                    #@assert test_quadratic(Λ_new, ρ, 0) == false
                end
            end

            if ζ_level_insertion == false
                break
            end
        end
        Λ[i] = Λ_new
    end

    return Λ
end


function recover_optimal_index_set{T}(Λ::PiecewiseQuadratic{T}, first_index=1)

    p, y, f = find_minimum(Λ)

    I = [first_index]

    while true
        push!(I, p.time_index)

        if !isdefined(p, :ancestor);
            break;
        end

        p = p.ancestor

        if p.time_index == -1
            break
        end
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


# TODO: Include mζ in the cost?!
"""
    I, Y, f = recover_solution(Λ::PiecewiseQuadratic{T}, ℓ, V_0N::QuadraticPolynomial, first_index=1)

"""
function recover_solution(Λ::PiecewiseQuadratic, ℓ, V_0N::QuadraticPolynomial, ζ=0.0)
    I, _, f_expected = recover_optimal_index_set(Λ, 1)
    Y, f = find_optimal_y_values(ℓ, V_0N::QuadraticPolynomial, I)

    f_regularized = f + ζ*(length(I)-1) # Include regularization cost

    !isapprox(f_regularized, f_expected, atol=1e-10) && warn("Recovered cost ($f_regularized) is not what was expected from value function ($f_expected). Solution might be incorrect.")

    return I, Y, f
end

"""
Evaluate the optimal cost (using least squares) for all
possible index sets with m segemets
"""
function brute_force_optimization(ℓ, V_0N::QuadraticPolynomial, m::Integer)
    cost_best = Inf

    I_best = []
    Y_best = []

    N = size(ℓ, 2)

    for I=IterTools.subsets(2:N-1, m-1)

        I = [1; I; N]
        P = zeros(m+1, m+1)
        q = zeros(m+1)
        r = 0

        # Add cost at right endpoint
        P[end,end] = V_0N.a
        q[end]     = V_0N.b
        r          = V_0N.c
        # Form quadratic cost function Y'*P*Y + q'*Y + r
        # corresponding to the y-values in the vector Y
        for j=1:m
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



"""
Given a piecewise quadratic V_Λ represented by Λ, a polynomial ρ, and a real number ζ
this function evaluates if
ρ(y) > V_Λ(y) + ζ  ∀ y
"""
# FIXME: Group polynomial and constant in tuple?
function poly_minus_constant_is_greater{T}(Λ::PiecewiseQuadratic{T}, ρ::QuadraticPolynomial{T}, ζ::Real)

    if Λ.next.left_endpoint == Inf
        return false
    end

    for λ_curr in Λ

        left_endpoint = λ_curr.left_endpoint
        right_endpoint = get_right_endpoint(λ_curr)

        Δa = ρ.a - λ_curr.p.a
        Δb = ρ.b - λ_curr.p.b
        Δc = (ρ.c - ζ) - λ_curr.p.c

        b2_minus_4ac =  Δb^2 - 4*Δa*Δc

        if Δa > 0 # ρ has greater curvature, i.e., ρ is smallest in the middle if intersect
            #println("Δa > 0")
            if b2_minus_4ac <= 0
                # Zero (or one) intersections, old quadratic is smallest, just step forward
                #println("No intersections, old quadratic is smallest, Δa > 0, breaking.")
                return true
            else

                # Compute the intersections
                term1 = -(Δb / 2 / Δa)
                term2 = sqrt(b2_minus_4ac) / 2 / Δa # Δa > 0
                root1, root2 = term1-term2, term1+term2

                # Check where the intersections are and act accordingly
                if root1 >= right_endpoint
                    continue
                elseif root2 <= left_endpoint
                    return true # There will be no more intersections since Δa > 0
                else
                    return false
                end
            end

        elseif Δa < 0 # ρ has lower curvature, i.e., ρ is smallest on the sides
            #println("Δa < 0")
            if b2_minus_4ac <= 0
                # Zero (or one) roots
                return false
            else
                # Compute the intersections
                term1 = -(Δb / 2 / Δa)
                term2 = sqrt(b2_minus_4ac) / 2 / Δa # < 0
                root1, root2 = term1+term2, term1-term2

                # Check where the intersections are and act accordingly
                if root1 <= left_endpoint && root2 >= right_endpoint
                    # One intersection on either side of the interval,
                    # old quadratic is smallest, just step forward
                    continue
                else
                    return false
                end
            end
        else # Δa == 0.0
            DEBUG2 && println("Δa == 0 : $ρ")
            if Δb == 0
                if Δc >= 0

                else
                    return false
                end
                continue
            end

            root = -Δc / Δb
            if Δb > 0
                if root < left_endpoint

                elseif root > right_endpoint
                    return false
                else
                    return false
                end
            else
                if root < left_endpoint
                    return false
                elseif root > right_endpoint

                else
                    DEBUG2 && println("Special case")
                    return false
                end
            end
        end

    end
    return true
end
