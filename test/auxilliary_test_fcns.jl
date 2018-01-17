# Some help functions to generate interesting test cases
circle_segment(N) = sin.(acos.(linspace(-1, 1, N)))
linear_trend(N) = (0:N-1)

# function get_transition_costs(g, t; tol=1e-3)
#     ℓ = TransitionCostContinuous{Float64}(g, t, tol=tol)
#     V_N = zero(QuadraticPolynomial{Float64})
#     return ℓ, V_N
# end
# function get_transition_costs(g::AbstractArray, t; tol=1e-3)
#     ℓ = compute_discrete_transition_costs(g, t)
#     V_N = QuadraticPolynomial(1.0, -2*g[t[end]], g[t[end]]^2)
#     return ℓ, V_N
# end



""" brute_force_multi(g, M, t=1:length(g); tol=1e-3)
    Finds the solution to the cardinality constrained problem for
    m=1:M segments using brute force with possible breakpoints on grid t.
    `t` has to be supplied if `g` is not `AbstractArray`
"""
function brute_force_multi(g, M, t=1:length(g); tol=1e-3)

    if length(g) <= M
        warn("Specified M too large. Setting M=length(g)-1.")
        M = length(g)-1
    end

    ℓ, V_N = DynamicApproximations.get_transition_costs(g, t, true, tol=tol)

    I_vec = Vector{Vector{Int64}}(M)
    Y_vec = Vector{Vector{Float64}}(M)
    f_vec = Vector{Float64}(M)

    for m=1:M
        (I_bf, Y_bf, f_bf) = brute_force_search(ℓ, V_N, m)
        I_vec[m] = t[I_bf]
        Y_vec[m] = Y_bf
        f_vec[m] = f_bf
    end

    return I_vec, Y_vec, f_vec
end
