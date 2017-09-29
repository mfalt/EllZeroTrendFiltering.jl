
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
    # T = Float64
    # # Find primitive functions to g, t*g, and g^2
    # # and evaluate them at the break points
    # #I_g = polyint(g).(t)
    # #I_g2 = polyint(g^2).(t)
    # #I_tg = polyint(Poly([0,1]) * g).(t)
    #
    # N = length(t)
    #
    # # Find primitive functions to g, t*g, and g^2 at the break points
    # I_g = zeros(N)
    # I_g2 = zeros(N)
    # I_tg = zeros(N)
    #
    # for i=2:N
    #     I_g[i] = I_g[i-1] + quadgk(g, t[i-1], t[i], reltol=tol)[1]
    #     I_g2[i] = I_g2[i-1] + quadgk(τ -> g(τ)^2, t[i-1], t[i], reltol=tol)[1]
    #     I_tg[i] = I_tg[i-1] + quadgk(τ -> τ*g(τ), t[i-1], t[i], reltol=tol)[1]
    # end
    #
    # ℓ = Array{QuadraticForm{T}}(N-1,N)
    #
    # for i=1:N-1
    #     for ip=i+1:N
    #
    #         P = (t[ip] - t[i]) * @SMatrix [1/3 1/6; 1/6 1/3]
    #
    #         q = -2* 1/(t[ip]-t[i]) *
    #         @SVector [-(I_tg[ip] - I_tg[i]) + t[ip]*(I_g[ip] - I_g[i]),
    #         (I_tg[ip] - I_tg[i]) - t[i]*(I_g[ip] - I_g[i])]
    #
    #         r =  I_g2[ip] - I_g2[i]
    #
    #         ℓ[i,ip] = QuadraticForm(P, q, r)
    #     end
    # end
    # TODO move/add tests to Lazy, maybe only i to i+1
    for i = 1:N-1, ip = i+1:N
        P, q, r = ℓ[i,ip].P, ℓ[i,ip].q, ℓ[i,ip].r
        minval = -q⋅(P\q)/4+r
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
