# Some help functions to generate interesting test cases
circle_segment(N) = sin.(acos.(linspace(-1, 1, N)))
linear_trend(N) = (0:N-1)

# Finds the solution to the cardinality constrained problem for
# m=1:M segments using brute force
function brute_force_multi(g, t, M; tol=1e-3)

    if length(g) <= M
        warn("Specified M too large. Setting M=length(g)-1.")
        M = length(g)-1
    end

    ℓ = []
    V_N = []

    if typeof(g) <: AbstractArray
        ℓ = compute_discrete_transition_costs(g)
        V_N = QuadraticPolynomial(1.0, -2*g[end], g[end]^2)
    else
        ℓ = TransitionCostContinuous{Float64}(g, t, tol=tol)
        V_N = zero(QuadraticPolynomial{Float64})
    end

    I_vec = Vector{Vector{Int64}}(M)
    Y_vec = Vector{Vector{Float64}}(M)
    f_vec = Vector{Float64}(M)

    for m=1:M
        (I_bf, Y_bf, f_bf) = brute_force_search(ℓ, V_N, m)
        I_vec[m] = I_bf
        Y_vec[m] = Y_bf
        f_vec[m] = f_bf
    end

    return I_vec, Y_vec, f_vec
end
