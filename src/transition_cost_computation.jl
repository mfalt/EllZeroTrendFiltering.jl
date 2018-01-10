
# TODO Maybe use big T for time indices, howabout mathcal{T}
"""
    ℓ = compute_transition_costs(g, t::AbstractArray)
Computes the transition costs `ℓ` given a
function `g` and a time sequence `t`
"""
function compute_transition_costs(g, t::AbstractVector, tol=1e-3)
    T = Float64
    ℓ_lazy = TransitionCostContinuous{T}(g, t)
    N = length(t)
    ℓ = Array{QuadraticForm{T}}(N-1,N)
    for i = 1:N-1
        for ip = i+1:N
            ℓ[i,ip] = ℓ_lazy[i,ip]
        end
    end
    # TODO move/add tests to Lazy, maybe only i to i+1
    for i = 1:N-1, ip = i+1:N
        P, q, r = ℓ[i,ip].P, ℓ[i,ip].q, ℓ[i,ip].r
        minval = -q⋅(P\q)/4+r
        #TODO add epsilon on RHS?
        minval < 0 && error("Transition cost is negative for some values in `compute_transition_costs`, try decreasing rtol.")
    end

    return ℓ
end

"""
    ℓ = compute_transition_costs(g::AbstractArray)
Computes the transition costs `ℓ` given a discrete
function `g`.
"""
function compute_discrete_transition_costs(g::AbstractArray; t=1:length(g))
    T = Float64
    ℓ_lazy = TransitionCostDiscrete{T}(g; t=t)
    N = length(t)
    ℓ = Array{QuadraticForm{T}}(N-1,N)
    for i = 1:N-1
        for ip = i+1:N
            ℓ[i,ip] = ℓ_lazy[i,ip]
        end
    end
    return ℓ
end
