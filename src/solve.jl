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

    Λ = Matrix{Union{PiecewiseQuadratic{Float64},Nothing}}(nothing, M, N)

    for i=1:N-1
        p = minimize_wrt_x2(l[i,N] + (V_N ∘ χ[i,N]))
        p.time_index = i
        p.ancestor = V_N
        Λ[1, i] = create_new_pwq(p)
    end

    μ = QuadraticPolynomial{T}()
    upper_bound_inner = Inf
    for m=2:M
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

    Λ = Vector{Union{PiecewiseQuadratic{T},Nothing}}(undef, N)

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




function recover_optimal_index_set_zero_ic(l::EllZeroTrendFiltering.AbstractTransitionCost{T}, Λ::Matrix{Union{PiecewiseQuadratic{T},Nothing}}, m::Integer) where T
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
	        end
	    end
	end
	I = recover_ancestors(λ_best.p)[1:end-1]
    I, cost_best
end


function recover_ancestors(p::QuadraticPolynomial)

    J = Vector{Int}(undef, 0)

    while true
        push!(J, p.time_index)

        if !isdefined(p, :ancestor);
            break;
        end

        p = p.ancestor
    end

    return J
end

function find_optimal_first_impulse(Λ::Vector{Union{PiecewiseQuadratic{T},Nothing}}, l) where T
    cost_best = typemax(T)
    p_best = QuadraticPolynomial{T}() # TODO: maybe could be QuadraticPolynomial instead

    for ip=2:length(Λ) # Or should it be N-m+1 or something else
        if Λ[ip] == nothing; break; end

        r = l[1,ip].r # Only handles zero initial conditions, arbitrary intiial conditions would be more messy
        for λ in Λ[ip]
            cost = r +  λ.p(0.0) # Incurred cost of doing nothing up to time ip, then impulse, x1 remains 0
            if cost_best >  cost
                cost_best = cost
                p_best = λ.p
            end
        end
    end

    return p_best, cost_best
end

function recover_optimal_index_set(Λ::Vector{Union{PiecewiseQuadratic{T},Nothing}}, l, χ, initial_conditions::Symbol=:zero) where T

    if initial_conditions == :free # free initial conditions
        p, y, f = find_minimum(Λ[1])

        # Find optimal initial conditions
        i = p.time_index
        ip = p.ancestor.time_index
        x0 = find_minimum(l[i, ip] + (p ∘ χ[i,ip]))[1]

        # Find optimal index set
        J = recover_ancestors(p)
        return J, f, x0, y
    elseif initial_conditions == :zero
        p_best, cost_best = find_optimal_first_impulse(Λ, l)

    	J = recover_ancestors(p_best)[1:end-1]
        return J, cost_best, [0.0; 0], 0.0
    else
        error("Unknown type of initial condition")
    end
end


function recover_optimal_index_set(Λ::PiecewiseQuadratic{T}) where T

    p, y, f = find_minimum(Λ)

    I = Vector{Int}(undef, 0)

    while true
        push!(I, p.time_index)

        if !isdefined(p, :ancestor); break; end
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

    if !isapprox(f_regularized, f_expected, atol=1e-8, rtol=1e-8)
        @warn "Recovered cost ($f_regularized) is not what was expected from value function ($f_expected). Solution might be incorrect."
    end

    if f_regularized < 0
        @warn "Computed cost ($f_regularized) < 0, if ≈ 0, this is probably due to numerical errors and nothing to worry about."
    end

    return I, Y, f
end




function post_process_pwl(l, V_N, J_vec::AbstractVector, f_expected_vec::AbstractVector, ζ=0.0)
    K = length(J_vec)
    Y_vec = Vector{Vector{Float64}}(undef, K)
    y_vec = Vector{Vector{Float64}}(undef, K)

    for k=1:K
        Y_vec, f = find_optimal_y_values(l, V_N::QuadraticPolynomial, J_vec[k])

        f_regularized = f + ζ*(length(I)-1) # Include regularization cost

        if !isapprox(f_regularized, f_expected_vec[k], atol=1e-8, rtol=1e-8)
            @warn "Recovered cost ($f_regularized) is not what was expected from value function ($f_expected). Solution might be incorrect."
        end

        if f_regularized < 0
            @warn "Computed cost ($f_regularized) < 0, if ≈ 0, this is probably due to numerical errors and nothing to worry about."
        end

        y_vec[k] = Y # FIXME, needs fixing
    end

    return Y_vec, y_vec
end
function post_process_lti(A::AbstractMatrix, C::AbstractMatrix, g, J_vec::AbstractVector; N_ic=0, initial_conditions=:zero)
    N = length(g)
    K = length(J_vec)
    U_vec = Vector{Vector{Float64}}(undef, K)
    y_vec = Vector{Vector{Float64}}(undef, K)
    x0_vec = Vector{Vector{Float64}}(undef, K)

    sys = ControlSystems.ss(Matrix{Float64}(A), [0.0; 1], Matrix{Float64}(C), 0, 1)
    X = generate_markov_matrix(sys, N; initial_conditions=initial_conditions)

    for k=1:K
        Xs = X[:, J_vec[k]]
        U_vec[k] = Xs\g
        y_vec[k] = Xs*U_vec[k]
    end

    return U_vec, y_vec # TODO: maybe also x0 out
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



function  fit_lti_constrained(A::AbstractMatrix, C::AbstractMatrix, g, M; tol=1e-3, precompute=false, initial_conditions=:zero)

    l, χ, V_N = compute_problem_data_lti(g, A, C)

    Λ = construct_value_fcn_constrained(l, χ, V_N, M, 1000)

    J_vec = Vector{Vector{Int}}(undef, M)
    f_vec = Vector{Float64}(undef, M)
    for m=1:M
        J_vec[m], f_vec[m] = recover_optimal_index_set(Λ[m, :], l, χ, initial_conditions)
    end
    J_vec = [J_vec[k] .- 1 for k=1:length(J_vec)]

    U_vec, y_vec = post_process_lti(A, C, g, J_vec, initial_conditions=initial_conditions)
    return J_vec, U_vec, f_vec, y_vec
end

function  ell0_constrained_dp(l::AbstractTransitionCost{T}, χ, V_N, M) where T
    Λ = construct_value_fcn_constrained(l, χ, V_N, M);

    I_vec = Vector{Vector{Int}}(undef, M)
    y0_vec = Vector{Vector{T}}(undef, M) # optimal first y-value
    f_vec = Vector{T}(undef, M)

    for m=1:M
        I_vec[m], y0_vec[m], f_vec[m] = recover_solution(Λ[m, 1], l, V_N)
    end

    return I_vec, y0_vec, f_vec
end
