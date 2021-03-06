# Some help functions to generate interesting test cases
circle_segment(N) = sin.(acos.(range(-1, stop=1, length=N)))
linear_trend(N) = (0:N-1)

""" brute_force_multi(g, M, t=1:length(g); tol=1e-3)
    Finds the solution to the cardinality constrained problem for
    m=1:M segments using brute force with possible breakpoints on grid t.
    `t` has to be supplied if `g` is not `AbstractArray`
"""
function brute_force_multi(g, M, t=1:length(g); tol=1e-3)

    if length(g) <= M
        @warn "Specified M too large. Setting M=length(g)-1."
        M = length(g)-1
    end

    ℓ, V_N = EllZeroTrendFiltering.get_transition_costs(g, t, true, tol=tol)

    I_vec = Vector{Vector{Int64}}(undef, M)
    Y_vec = Vector{Vector{Float64}}(undef, M)
    f_vec = Vector{Float64}(undef, M)

    for m=1:M
        (I_bf, Y_bf, f_bf) = brute_force_search(ℓ, V_N, m)
        I_vec[m] = t[I_bf]
        Y_vec[m] = Y_bf
        f_vec[m] = f_bf
    end

    return I_vec, Y_vec, f_vec
end
