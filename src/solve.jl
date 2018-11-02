
"""
    add_quadratic!(Λ::PiecewiseQuadratic{T}, μ::QuadraticPolynomial{T}) where {T}
Inserts a quadratic polynomial μ into the linked list `Λ`, which represents a piecewise
quadratic polynomial, so that the new `Λ` satisfies
`Λ(y) := min{ Λ(y) ,μ(y) } ∀ y`
"""
function add_quadratic!(Λ::PiecewiseQuadratic{T}, μ::QuadraticPolynomial{T}) where T

    if Λ.next.left_endpoint == Inf # I.e. the piecewise quadratic object is empty, perhaps better to add dummy polynomial
        μ.has_been_used = true
        insert(Λ, μ, -Inf) # FIXME: Is specific insert function needed?
        return
    end

    λ_prev = Λ # Points to the list head which is just a node with NaN data
    λ_curr = Λ.next # The leftmost segment of the list

    while λ_curr.left_endpoint != Inf #???
        #global counter2 += 1
        DEBUG && println(Λ)

        left_endpoint = λ_curr.left_endpoint
        right_endpoint = get_right_endpoint(λ_curr)


        Δa = μ.a - λ_curr.p.a
        Δb = μ.b - λ_curr.p.b
        Δc = μ.c - λ_curr.p.c

        b2_minus_4ac =  Δb^2 - 4*Δa*Δc

        if Δa > 0 # μ has greater curvature, i.e., μ is smallest in the middle if intersect
            if b2_minus_4ac <= ACCURACY
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
                    DEBUG && println("Two intersections to the right: old quadratic is smallest")
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                elseif root2 <= left_endpoint
                    DEBUG && println("Two intersections to the left: old quadratic is smallest")
                    break # There will be no more intersections since Δa > 0
                elseif root1 <= left_endpoint && root2 >= right_endpoint
                    DEBUG && println("One intersections on either side, new quadratic is smallest")
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
                elseif root1 > left_endpoint && root2 < right_endpoint
                    DEBUG && println("Two intersections: old-new-old")
                    λ_prev, λ_curr = update_segment_old_new_old(λ_curr, μ, root1, root2)
                    break # There will be no more intersections since Δa > 0
                elseif root1 > left_endpoint
                    DEBUG && println("Root 1 within the interval: old-new")
                    λ_prev, λ_curr = update_segment_old_new(λ_curr, μ, root1)
                elseif root2 < right_endpoint
                    DEBUG && println("Root 2 within the interval: new-old")
                    λ_prev, λ_curr = update_segment_new_old(λ_prev, λ_curr, μ, root2)
                    break # There will be no more intersections since Δa > 0
                else
                    error("Shouldn't end up here")
                end
            end

        elseif Δa < 0 # μ has lower curvature, i.e., μ is smallest on the sides
            if b2_minus_4ac <= ACCURACY
                # Zero (or one) roots
                λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
            else
                # Compute the intersections
                term1 = -(Δb / 2 / Δa)
                term2 = sqrt(b2_minus_4ac) / 2 / Δa # < 0
                root1, root2 = term1+term2, term1-term2
                DEBUG && println("Δa < 0   root1:", root1, "   root2:", root2)

                # Check where the intersections are and act accordingly
                if root1 >= right_endpoint || root2 <= left_endpoint
                    # No intersections, μ is smallest
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
                elseif root1 <= left_endpoint && root2 >= right_endpoint
                    # No intersections, old quadratic is smallest, just step forward
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                elseif root1 > left_endpoint && root2 < right_endpoint
                    # Two intersections within the interval
                    λ_prev, λ_curr = update_segment_new_old_new(λ_prev, λ_curr, μ, root1, root2)
                elseif root1 > left_endpoint
                    λ_prev, λ_curr = update_segment_new_old(λ_prev, λ_curr, μ, root1)
                elseif root2 < right_endpoint
                    λ_prev, λ_curr = update_segment_old_new(λ_curr, μ, root2)
                else
                    error("Shouldn't end up here")
                end
            end
        else # Δa == 0.0
            DEBUG && println("Δa == 0")
            DEBUG2 && println("Δa == 0 : $μ")
            if Δb == 0
                if Δc >= 0
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                else
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
                end
                continue
            end

            root = -Δc / Δb
            if Δb > 0
                if root < left_endpoint
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                elseif root > right_endpoint
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
                else
                    λ_prev, λ_curr = update_segment_new_old(λ_prev, λ_curr, μ, root)
                end
            else
                if root < left_endpoint
                    λ_prev, λ_curr = update_segment_new(λ_prev, λ_curr, μ)
                elseif root > right_endpoint
                    λ_prev, λ_curr = update_segment_do_nothing(λ_curr)
                else
                    DEBUG2 && println("Special case")
                    λ_prev, λ_curr = update_segment_old_new(λ_curr, μ, root)
                end
            end
        end

    end
    return
end

@inline function update_segment_do_nothing(λ_curr)
    return λ_curr, λ_curr.next
end

@inline function update_segment_new(λ_prev, λ_curr, μ)
    μ.has_been_used = true # TODO: Could be moved to the second clause
    if λ_prev.p === μ
        λ_prev.next = λ_curr.next
        v1, v2 = λ_prev, λ_curr.next
    else
        λ_curr.p = μ
        v1, v2 = λ_curr, λ_curr.next
    end
    return v1, v2
end

@inline function update_segment_old_new(λ_curr, μ, root)
    # ... λ_curr | λ_new ...
    μ.has_been_used = true
    λ_new =  PiecewiseQuadratic(μ, root, λ_curr.next)
    λ_curr.next = λ_new
    return λ_new, λ_new.next
end

@inline function update_segment_new_old(λ_prev, λ_curr, μ, root)
    # ... λ_prev | (new segment) | λ_curr ...
    if λ_prev.p === μ
        λ_curr.left_endpoint = root
    else
        μ.has_been_used = true
        λ_prev.next = PiecewiseQuadratic(μ, λ_curr.left_endpoint, λ_curr)
        λ_curr.left_endpoint = root
    end
    return λ_curr, λ_curr.next
end

@inline function update_segment_new_old_new(λ_prev, λ_curr, μ, root1, root2)
    # Insert new segments with μ before and after λ_curr
    # ... λ_prev : (new segment) | λ_curr | (new segment) ...
    # ( μ.has_been_used = true is set inside the called funcitons )
    update_segment_new_old(λ_prev, λ_curr, μ, root1)
    return update_segment_old_new(λ_curr, μ, root2)
end

@inline function update_segment_old_new_old(λ_curr, μ, root1, root2)
    # ... λ_curr | λ_1  | λ_2 ...
    # The second new piece contains a copy of the polynomial in λ_curr

    μ.has_been_used = true
    λ_2 =  PiecewiseQuadratic(λ_curr.p, root2, λ_curr.next)
    λ_1 =  PiecewiseQuadratic(μ, root1, λ_2)
    λ_curr.next = λ_1
    return λ_2, λ_2.next
end

"""
    μ = minimize_wrt_x2(qf::QuadraticForm{T}, p::QuadraticPolynomial{T}, μ=QuadraticPolynomial{T}()) where {T}
Takes a quadratic form in `[x₁; x₂]` and a polynomial in `x₂`
and returns the minimum of the sum wrt to `x₂`,
i.e. μ(x₁) = min_x₂{ qf(x₁,x₂) +  p(x₂) }`

The input `μ` can be pre-allocated on input and will then be changed.
"""
@inline function minimize_wrt_x2_old(qf::QuadraticForm{T}, p::QuadraticPolynomial{T}, μ=QuadraticPolynomial()) where T
    P = qf.P
    q = qf.q
    r = qf.r

    P22_new = P[2,2] + p.a

    if P22_new > 0
        μ.a = P[1,1] - P[1,2]^2 / P22_new
        μ.b = q[1] - P[1,2]*(q[2] + p.b) / P22_new
        μ.c = (r + p.c) - (q[2] + p.b)^2 / P22_new / 4
    elseif P22_new == 0
        μ.a = P[1,1]
        μ.b = q[1]
        μ.c = r + p.c
    else
        error("Piecewise quadratic cost should be positive (semi-)definite. Please submit a bug report.")
    end
    return μ
end


@inline function minimize_wrt_x2_alt(qf::QuadraticForm{T}, χ::SVector{2,T}, p::QuadraticPolynomial{T}, μ=QuadraticPolynomial{T}()) where T
    P = qf.P + p.a*χ*χ'
    q = qf.q + p.b*χ
    r = qf.r + p.c

    if P[2,2] > 0
        μ.a = P[1,1] - P[1,2]^2 / P[2,2]
        μ.b = q[1] - P[1,2]*q[2] / P[2,2]
        μ.c = r - q[2]^2 / P[2,2] / 4
    elseif P[2,2] == 0
        μ.a = P[1,1]
        μ.b = q[1]
        μ.c = r
    else
        error("Piecewise quadratic cost should be positive (semi-)definite. Please submit a bug report.")
    end
    return μ
end


@inline function minimize_wrt_x2(qf::QuadraticForm{T}, μ=QuadraticPolynomial{T}()) where T
    P = qf.P
    q = qf.q
    r = qf.r

    if P[2,2] > 0
        μ.a = P[1,1] - P[1,2]^2 / P[2,2]
        μ.b = q[1] - P[1,2]*q[2] / P[2,2]
        μ.c = r - q[2]^2 / P[2,2] / 4
    elseif P[2,2] == 0
        μ.a = P[1,1]
        μ.b = q[1]
        μ.c = r
    else
        error("Piecewise quadratic cost should be positive (semi-)definite. Please submit a bug report.")
    end
    return μ
end





global counter1
global counter2

"""
Construct the value function that corresponds to the ell_0 constrained problem

minimize   ∑        l[i, i'](y, y')  +  V_N(χ[i,i'] y)  +  ζ⋅
         i,i' ∈ I

subject to  card(I) ≦ M

where l[j,k] are positive-definite quadratic forms.

Note that the value functions corresponding to card(I) ≦ m, m < M are also returned.
"""
function construct_value_fcn_constrained(l::AbstractTransitionCost{T}, χ::AbstractMatrix, V_N::QuadraticPolynomial{T}, M::Integer, upper_bound=Inf) where T
    #global counter1
    #global counter2
    #counter1 = 0
    #counter2 = 0

    N = size(l, 2)
    @assert M <= N-1 "Cannot have more segments than N-1."

    Λ = Array{PiecewiseQuadratic{T}}(undef, M, N)

    for i=1:N-1
        p = minimize_wrt_x2(l[i,N] + (V_N ∘ χ[i,N]))
        p.time_index = N
        Λ[1, i] = create_new_pwq(p)
    end

    μ = QuadraticPolynomial{T}()
    upper_bound_inner = Inf
    global times
    for m=2:M
        #println("m: $m")
        for i=1:N-m
            Λ_new = create_new_pwq(T)
            if OPTIMIZE
                if min(upper_bound, upper_bound_inner) < Inf
                    add_quadratic!(Λ_new, QuadraticPolynomial{T}(0.0, 0.0, min(upper_bound,upper_bound_inner)))
                end
            end

            for ip=i+1:N-m+1
                DEBUG && println("(m:$m, i:$i, ip:$ip)")

                for λ in Λ[m-1, ip]

                    minimize_wrt_x2(l[i,ip] + (λ.p ∘ χ[i,ip]), μ)

                    DEBUG && println("Obtained μ = $μ")

                    # if OPTIMIZE
                    #     μmin = unsafe_minimum(μ)
                    #     if (μmin > upper_bound || μmin > upper_bound_inner)
                    #         DEBUG && println("Breaking due to that $μmin > max($upper_bound, $upper_bound_inner)")
                    #         continue
                    #     end
                    # end

                    DEBUG && println("Inserting...")

                    add_quadratic!(Λ_new, μ)

                    if μ.has_been_used == true
                        μ.time_index = ip
                        μ.ancestor = λ.p
                        μ = QuadraticPolynomial{T}()
                        μ.has_been_used = false
                    end
                end
            end
            #remove_over(Λ_new, min(upper_bound_inner,upper_bound))
            Λ[m, i] = Λ_new
            if OPTIMIZE
                if i == 1
                    upper_bound_inner = find_minimum_value(Λ[m,1])
                    upper_bound_inner += sqrt(eps())
                end
            end
        end
    end
    return Λ
end



"""
Construct the value function that corresponds to the regularization problem

minimize   ∑        l[i, i'](y, y')  +  V_N(χ[i,i'] y)  +  ζ⋅card(I)
         i,i' ∈ I
where l[j,k] are positive-definite quadratic forms.
"""
function construct_value_fcn_regularized(l::AbstractTransitionCost{T}, χ::AbstractMatrix, V_N::QuadraticPolynomial{T}, ζ::T) where T
    N = size(l, 2)

    Λ = Vector{PiecewiseQuadratic{T}}(undef, N)

    V_N = deepcopy(V_N)
    V_N.time_index = -1
    Λ[N] = create_new_pwq(V_N)

    μ = QuadraticPolynomial{T}()

    for i=N-1:-1:1
        Λ_new = create_new_pwq(T)
        for ip=i+1:N


            ζ_level_insertion = false
            for λ in Λ[ip]

                minimize_wrt_x2(l[i,ip] + (λ.p ∘ χ[i,ip]), μ)

                μ.c += ζ # add cost for break point


                add_quadratic!(Λ_new, μ)

                if μ.has_been_used
                    μ.time_index = ip
                    μ.ancestor = λ.p
                    μ = QuadraticPolynomial{T}()
                    μ.has_been_used = false

                    ζ_level_insertion = true
                else
                    if ζ_level_insertion == false
                        if !poly_minus_constant_is_greater(Λ_new, μ, ζ)
                            ζ_level_insertion = true
                        end
                    end
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




function recover_optimal_index_set_zero_ic(l::EllZeroTrendFiltering.AbstractTransitionCost{T}, Λ::Matrix{PiecewiseQuadratic{T}}, m::Integer) where T
    N = size(l, 2)
	cost_best = 10000000
	λ_best = -1
    ip_best = -1
	# For the constrained case
	for ip=2:N-m # Or should it be N-m+1 or something else
		r = l[1,ip].r # Only handles zero initial conditions, arbitrary intiial conditions would be more messy
	    for λ in Λ[m, ip]
	        if cost_best > r +  λ.p(0.0)
	            cost_best = r + λ.p(0.0)
                ip_best = ip
	            λ_best = λ
                println("$(λ_best.p) : $cost_best")
	        end
	    end
	end
	I = [ip_best; recover_ancestors(λ_best.p)[1:end-1]]
    I, cost_best
end

function recover_ancestors(p::QuadraticPolynomial)
	I = Vector{Int}(undef, 0)

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

	return I
end


function recover_optimal_index_set(Λ::PiecewiseQuadratic{T}, first_index=1) where T

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
    Y, f = find_optimal_y_values(l, V_N, I)
Given transition costs `l`, cost of right endpoint `V_N`, and breakpoint indicies `I`
the optimal y-values `Y` and the optimal cost `f` are computed.
"""
function find_optimal_y_values(l, V_N::QuadraticPolynomial, I)
    m = length(I) - 1

    P = zeros(m+1, m+1)
    q = zeros(m+1)

    # Add cost for the right endpoint
    P[m+1, m+1] = V_N.a
    q[m+1] = V_N.b
    r = V_N.c

    # Form quadratic cost function Y'*P*Y + q'*Y + r
    # corresponding to the y-values in the vector Y
    for j=1:m
        P[j:j+1,j:j+1] .+= l[I[j], I[j+1]].P
        q[j:j+1] .+= l[I[j], I[j+1]].q
        r += l[I[j], I[j+1]].r
    end

    # find the optimal Y-vector, and compute the correspinding error
    Y = -(P \ q) / 2
    f = Y' * P * Y + q' * Y + r
    return Y, f
end


# TODO: Include mζ in the cost?!
"""
    I, Y, f = recover_solution(Λ::PiecewiseQuadratic{T}, l, V_N::QuadraticPolynomial, first_index=1)

"""
function recover_solution(Λ::PiecewiseQuadratic, l, V_N::QuadraticPolynomial, ζ=0.0)
    I, _, f_expected = recover_optimal_index_set(Λ, 1)
    Y, f = find_optimal_y_values(l, V_N::QuadraticPolynomial, I)

    f_regularized = f + ζ*(length(I)-1) # Include regularization cost

    if !isapprox(f_regularized, f_expected, atol=1e-10)
        @warn "Recovered cost ($f_regularized) is not what was expected from value function ($f_expected). Solution might be incorrect."
    end

    if f_regularized < 0
        @warn "Computed cost ($f_regularized) < 0, if ≈ 0, this is probably due to numerical errors and nothing to worry about."
    end

    return I, Y, f
end




"""
Given a piecewise quadratic V_Λ represented by Λ, a polynomial μ, and a real number ζ
this function evaluates if
μ(y) > V_Λ(y) + ζ  ∀ y
""" # FIXME: Group polynomial and constant in tuple?
function poly_minus_constant_is_greater(Λ::PiecewiseQuadratic{T}, μ::QuadraticPolynomial{T}, ζ::Real) where T

    if Λ.next.left_endpoint == Inf
        return false
    end

    for λ_curr in Λ

        left_endpoint = λ_curr.left_endpoint
        right_endpoint = get_right_endpoint(λ_curr)

        Δa = μ.a - λ_curr.p.a
        Δb = μ.b - λ_curr.p.b
        Δc = (μ.c - ζ) - λ_curr.p.c

        b2_minus_4ac =  Δb^2 - 4*Δa*Δc

        if Δa > 0 # μ has greater curvature, i.e., μ is smallest in the middle if intersect
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

        elseif Δa < 0 # μ has lower curvature, i.e., μ is smallest on the sides
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
                    # No intersection on either side of the interval, and μ is smallest
                    continue
                else
                    return false
                end
            end
        else # Δa == 0.0
            DEBUG2 && println("Δa == 0 : $μ")
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



## Public interface
"""
    I, Y, v = fit_pwl_regularized(g::AbstractArray, ζ; t=1:length(g), precompute=false)
    I, Y, v = fit_pwl_regularized(g, t, ζ, tol=1e-3; precompute=false)
Approximate `g[k]` (or `g(t)`) with a continuous piecewise linear function `f` according to
    v = min_f `||f-g||₂^2 + ζ⋅(length(I)-2)`
where the norm is `sum((f[k]-g[k])²)` (or integral over `(f(t)-g(t))²` from `t[1]` to `t[end]`).
Returns:
    `I`: Vector of length `M`
    `Y`: Vector of length `M`
    `v`: Float64
such that the optimal function has breakpoints in `I` (or `t[I]`) and satisfies
`f(I) .== Y` (or `f(t[I]) .== Y`).
Kwargs:
`t` is optional parameter in the discrete case, restricting the set of possible gridpoints,
i.e. so that `f[t[I]] .== Y`.
`precompute` = true, means that the internal transition costs `l[i,j]` will be calculated when needed.
`tol` specifies the relative tolerance sent to `quadg` kused when calculating the integrals (continuous case).
"""
function  fit_pwl_regularized(g::AbstractArray, ζ; precompute=false, t=1:length(g))
    l, χ, V_N = compute_problem_data_pwl(g, t; precompute=precompute)
    ell0_regularized_dp(l, χ, V_N, ζ)
end

# Continuous function
function  fit_pwl_regularized(g, t, ζ; tol=1e-3, precompute=false)
    l, χ, V_N = compute_problem_data_pwl(g, t; precompute=precompute, tol=tol)
    ell0_regularized_dp(l, χ, V_N, ζ)
end

function  fit_lti_output_regularized(sys::ControlSystems.StateSpace, g::AbstractArray, ζ; t=1:length(g), precompute=false)
    l, χ, V_N = compute_problem_data_lti(g, sys.A, sys.C) # TODO: fix arbitrary times
    ell0_regularized_dp(l, χ, V_N, ζ)
end
function  fit_lti_output_regularized(sys::ControlSystems.StateSpace, g, t, ζ, tol=1e-3; precompute=false)
    error("Not handled yet")
end

function  ell0_regularized_dp(l, χ, V_N, ζ)
    Λ_reg = construct_value_fcn_regularized(l, χ, V_N, ζ)
    #Get solution that starts at first index
    I, Y, f = recover_solution(Λ_reg[1], l, V_N, ζ)
    return I, Y, f
end


"""
    I, Y, v = fit_pwl_constrained(g::AbstractArray, M; t=1:length(g), precompute=false)
    I, Y, v = fit_pwl_constrained(g, t, M, tol=1e-3; precompute=false)
Approximate `g[k]` (or `g(t)`) with a continuous piecewise linear function `f` according to
    v = min_f `||f-g||₂^2`
    s.t. length(I)-2 = M
where the norm is `sum((f[k]-g[k])²)` (or integral over `(f(t)-g(t))²` from `t[1]` to `t[end]`).
Returns:
    `I`: Vector of length `M`
    `Y`: Vector of length `M`
    `v`: Float64
such that the optimal function has breakpoints in `I` (or `t[I]`) and satisfies
`f(I) .== Y` (or `f(t[I]) .== Y`).
Kwargs:
`t` is optional parameter in the discrete case, restricting the set of possible gridpoints,
i.e. so that `f[t[I]] .== Y`.
`precompute` = true, means that the internal transition costs `l[i,j]` will be calculated when needed.
`tol` specifies the relative tolerance sent to `quadg` kused when calculating the integrals (continuous case).
"""
function  fit_pwl_constrained(g::AbstractArray, M; t=1:length(g), precompute=true)
    l, χ, V_N = compute_problem_data_pwl(g, t; precompute=precompute)
    ell0_constrained_dp(l, χ, V_N, M)
end

# Continuous function
function  fit_pwl_constrained(g, t, M; tol=1e-3, precompute=false)
    l, χ, V_N = compute_problem_data_pwl(g, t; precompute=precompute, tol=tol)
    ell0_constrained_dp(l, χ, V_N, M)
end

function  ell0_constrained_dp(l::AbstractTransitionCost{T}, χ, V_N, M) where T
    Λ = construct_value_fcn_constrained(l, χ, V_N, M);

    Ivec = Vector{Vector{Int}}(undef, M)
    Yvec = Vector{Vector{T}}(undef, M)
    fvec = Vector{T}(undef, M)

    for m=1:M
        Ivec[m], Yvec[m], fvec[m] = recover_solution(Λ[m, 1], l, V_N)
    end

    return Ivec, Yvec, fvec
end
