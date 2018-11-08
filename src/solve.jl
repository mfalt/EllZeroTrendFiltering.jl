"""
Construct the value function that corresponds to the ell_0 constrained problem

minimize   ∑        l[i, i'](y, y')  +  V_N(χ[i,i'] y)  +  ζ
         i,i' ∈ I

subject to  card(I) ≦ M

where l[j,k] are positive-definite quadratic forms.

Note that the value functions corresponding to card(I) ≦ m, m < M are also returned.
"""
function construct_value_fcn_constrained(l::AbstractTransitionCost{T}, χ::AbstractMatrix, V_N::QuadraticPolynomial{T}, M::Integer, upper_bound=Inf) where T
    N = size(l, 2)

    V_N = deepcopy(V_N)
    V_N.time_index = N

    @assert M <= N-1 "Cannot have more segments than N-1."

    Λ = Array{PiecewiseQuadratic{T}}(undef, M, N)

    for i=1:N-1
        p = minimize_wrt_x2(l[i,N] + (V_N ∘ χ[i,N]))
        p.time_index = i
        p.ancestor = V_N
        Λ[1, i] = create_new_pwq(p)
    end

    μ = QuadraticPolynomial{T}()
    upper_bound_inner = Inf
    for m=2:M
        #println("m: $m")
        for i=1:N-m
            Λ_new = create_new_pwq(T)
            if OPTIMIZE
                if min(upper_bound, upper_bound_inner) < Inf
                    insert_quadratic!(Λ_new, QuadraticPolynomial{T}(0.0, 0.0, min(upper_bound,upper_bound_inner)))
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

                    insert_quadratic!(Λ_new, μ)

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
Construct the value function that corresponds to the regularization problem

minimize   ∑        l[i, i'](y, y')  +  V_N(χ[i,i'] y)  +  ζ⋅card(I)
         i,i' ∈ I
where l[j,k] are positive-definite quadratic forms.
"""
function construct_value_fcn_regularized(l::AbstractTransitionCost{T}, χ::AbstractMatrix, V_N::QuadraticPolynomial{T}, ζ::T) where T
    N = size(l, 2)

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

                minimize_wrt_x2(l[i,ip] + (λ.p ∘ χ[i,ip]), μ)

                μ.c += ζ # add cost for break point


                insert_quadratic!(Λ_new, μ)

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
	I = recover_ancestors(λ_best.p)[1:end-1]
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
	end

	return I
end


function recover_optimal_index_set(Λ::PiecewiseQuadratic{T}) where T

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
    I, _, f_expected = recover_optimal_index_set(Λ)
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
