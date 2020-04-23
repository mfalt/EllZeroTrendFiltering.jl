
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
@inline function minimize_wrt_x2(qf::QuadraticForm{T}, p::QuadraticPolynomial{T}, μ=QuadraticPolynomial{T}()) where T
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




global counter1
global counter2

"""
    Λ = pwq_dp_constrained(ℓ::AbstractTransitionCost{T}, V_N::QuadraticPolynomial{T}, M::Integer, upper_bound=Inf) where {T}

Given the transition costs `ℓ[i,j](y_i,y_j)` and the cost at the endpoint `V_N(y_N)` find all solutions `f` with up to `M` segments for the problem

V_i^m = minimize_f^M [ Σ_{k=1}^i { ℓ[k,k+1](f(k),f(k+1)) } + V_N(f(N)) ]
s.t          f(k) being continuous piecewise linear with `m` segements.

i.e. `Λ[m,i]` contains the best (in `ℓ` cost) continuous piecewise linear function `f` with up to `M` segments over the interval `i` to `N`
"""
function pwq_dp_constrained(ℓ::AbstractTransitionCost{T}, V_N::QuadraticPolynomial{T}, M::Integer, upper_bound=Inf) where T
    #global counter1
    #global counter2
    #counter1 = 0
    #counter2 = 0
    N = size(ℓ, 2)
    @assert M <= N-1 "Cannot have more segments than N-1."

    V_N = deepcopy(V_N)
    V_N.time_index = N

    Λ = Array{PiecewiseQuadratic{T}}(undef, M, N)

    for i=1:N-1
        p = minimize_wrt_x2(ℓ[i, N], V_N)
        p.time_index = i
        p.ancestor = V_N
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


                    minimize_wrt_x2(ℓ[i,ip], λ.p, μ)

                    DEBUG && println("Obtained μ = $μ")

                    if OPTIMIZE
                        μmin = unsafe_minimum(μ)
                        if (μmin > upper_bound || μmin > upper_bound_inner)
                            DEBUG && println("Breaking due to that $μmin > max($upper_bound, $upper_bound_inner)")
                            continue
                        end
                    end

                    DEBUG && println("Inserting...")

                    add_quadratic!(Λ_new, μ)

                    if μ.has_been_used == true
                        μ.time_index = i
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
TODO Update docstring:

Finds the set I=(i_1, ..., i_M) that is the solution to the regularization problem

minimize ∑ ℓ_(i, i+1)(y, y+1)  +  V_N(y_M)  +  ζ⋅card(I)

where ℓ are positive-definite quadratic forms.
"""
function pwq_dp_regularized(ℓ::AbstractTransitionCost{T}, V_N::QuadraticPolynomial{T}, ζ::T) where T
    N = size(ℓ, 2)

    Λ = Vector{PiecewiseQuadratic{T}}(undef, N)

    V_N = deepcopy(V_N)
    V_N.time_index = N
    Λ[N] = create_new_pwq(V_N)

    μ = QuadraticPolynomial{T}()

    for i=N-1:-1:1
        Λ_new = create_new_pwq(T)
        for ip=i+1:N


            ζ_level_insertion = false
            for λ in Λ[ip]
                #counter1 += 1

                minimize_wrt_x2(ℓ[i,ip], λ.p, μ)
                μ.c += ζ # add cost for break point


                add_quadratic!(Λ_new, μ)

                if μ.has_been_used
                    μ.time_index = i
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


function recover_optimal_index_set(Λ::PiecewiseQuadratic{T}, first_index=1) where T

    p, y, f = find_minimum(Λ)

    I = Vector{Int}(undef, 0)

    while true
        push!(I, p.time_index)

        if !isdefined(p, :ancestor);
            break;
        end

        p = p.ancestor
    end

    return I, y, f
end

"""
    Y, f = find_optimal_y_values(ℓ, V_N, I)
Given transition costs `ℓ`, cost of right endpoint `V_N`, and breakpoint indicies `I`
the optimal y-values `Y` and the optimal cost `f` are computed.
"""
function find_optimal_y_values(ℓ, V_N::QuadraticPolynomial, I)
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
    I, Y, f = recover_solution(Λ::PiecewiseQuadratic{T}, ℓ, V_N::QuadraticPolynomial, ζ=0.0)

"""
function recover_solution(Λ::PiecewiseQuadratic, ℓ, V_N::QuadraticPolynomial, ζ=0.0)
    I, _, f_expected = recover_optimal_index_set(Λ, 1)
    Y, f = find_optimal_y_values(ℓ, V_N::QuadraticPolynomial, I)

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

# Continuous case
""" `get_transition_costs(g, t, lazy; tol=1e-3)`
    Return `(ℓ, cost_last)`: the transition cost and end-cost.
    Slightly unsafe since `tol` is ignored (not needed) on discrete problems
"""
function get_transition_costs(g, t, lazy; tol=1e-3)
    T = Float64
    ℓ = lazy ? TransitionCostContinuous{T}(g, t, tol) :
               compute_transition_costs(g, t, tol)
    # Continouous case, no cost at endpoint
    cost_last = zero(QuadraticPolynomial{Float64})
    return ℓ, cost_last
end

# Discrete case, tol is not used here, but in signature to enable dispatch
function get_transition_costs(g::AbstractArray, t, lazy; tol=1e-3)
    T = promote_type(eltype(g), Float64)
    ℓ = lazy ? TransitionCostDiscrete{T}(g, t) :
               compute_discrete_transition_costs(g, t)
    if t[1] != 1 || t[end] != length(g)
        @warn "In fit_pwl_constrained: The supplied grid t only covers the range ($(t[1]),$(t[end])) while the range of indices for g is (1,$(length(g))). No costs will be considered outside the range of t."
    end
    # Discrete case, so cost at endpoint is quadratic
    cost_last = QuadraticPolynomial(1.0, -2*g[t[end]], g[t[end]]^2)
    return ℓ, cost_last
end


## Public interface
"""
    I, Y, v = fit_pwl_regularized(g::AbstractArray, ζ; t=1:length(g), lazy=true)
    I, Y, v = fit_pwl_regularized(g, t, ζ, tol=1e-3; lazy=true)
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
`lazy` = true, means that the internal transition costs `ℓ[i,j]` will be calculated when needed.
`tol` specifies the relative tolerance sent to `quadg` kused when calculating the integrals (continuous case).
"""
function  fit_pwl_regularized(g::AbstractArray, ζ; t=1:length(g), lazy=true)
    ℓ, cost_last = get_transition_costs(g, t, lazy)
    fit_pwl_regularized_internal(ℓ, cost_last, ζ)
end

# Continuous function
function  fit_pwl_regularized(g, t, ζ, tol=1e-3; lazy=true)
    ℓ, cost_last = get_transition_costs(g, t, lazy, tol=tol)
    fit_pwl_regularized_internal(ℓ, cost_last, ζ)
end

function  fit_pwl_regularized_internal(ℓ, cost_last, ζ)
    Λ_reg = pwq_dp_regularized(ℓ, cost_last, ζ)
    #Get solution that starts at first index
    I, Y, f = recover_solution(Λ_reg[1], ℓ, cost_last, ζ)
    return I, Y, f
end


"""
    I, Y, v = fit_pwl_constrained(g::AbstractArray, M; t=1:length(g), lazy=true)
    I, Y, v = fit_pwl_constrained(g, t, M, tol=1e-3; lazy=true)
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
`lazy` = true, means that the internal transition costs `ℓ[i,j]` will be calculated when needed.
`tol` specifies the relative tolerance sent to `quadg` kused when calculating the integrals (continuous case).
"""
function  fit_pwl_constrained(g::AbstractArray, M; t=1:length(g), lazy=false)
    ℓ, cost_last = get_transition_costs(g, t, lazy)
    fit_pwl_constrained_internal(ℓ, cost_last, M)
end

# Continuous function
function  fit_pwl_constrained(g, t, M, tol=1e-3; lazy=false)
    ℓ, cost_last = get_transition_costs(g, t, lazy, tol=tol)
    fit_pwl_constrained_internal(ℓ, cost_last, M)
end

function  fit_pwl_constrained_internal(ℓ::AbstractTransitionCost{T}, cost_last, M) where T
    Λ = pwq_dp_constrained(ℓ, cost_last, M);

    Ivec = Vector{Vector{Int}}(undef, M)
    Yvec = Vector{Vector{T}}(undef, M)
    fvec = Vector{T}(undef, M)

    for m=1:M
        Ivec[m], Yvec[m], fvec[m] = recover_solution(Λ[m, 1], ℓ, cost_last)
    end

    return Ivec, Yvec, fvec
end
