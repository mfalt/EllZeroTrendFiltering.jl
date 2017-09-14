
# TODO Maybe use big T for time indices, howabout mathcal{T}
"""
    ℓ = compute_transition_costs(g, t::AbstractArray)
Computes the transition costs `ℓ` given a
function `g` and a time sequence `t`
"""
function compute_transition_costs(g, t::AbstractArray)
    T = Float64
    # Find primitive functions to g, t*g, and g^2
    # and evaluate them at the break points
    #I_g = polyint(g).(t)
    #I_g2 = polyint(g^2).(t)
    #I_tg = polyint(Poly([0,1]) * g).(t)

    N = length(t)

    # Find primitive functions to g, t*g, and g^2 at the break points
    I_g = zeros(size(t))
    I_g2 = zeros(size(t))
    I_tg = zeros(size(t))

    for i=2:N
        const tol = 1e-3
        I_g[i] = I_g[i-1] + quadgk(g, t[i-1], t[i], reltol=tol)[1]
        I_g2[i] = I_g2[i-1] + quadgk(t -> g(t)^2, t[i-1], t[i], reltol=tol)[1]
        I_tg[i] = I_tg[i-1] + quadgk(t -> t*g(t), t[i-1], t[i], reltol=tol)[1]
    end

    ℓ = Array{QuadraticForm{T}}(N-1,N)

    for i=1:N-1
        for ip=i+1:N

            P = (t[ip] - t[i]) * @SMatrix [1/3 1/6; 1/6 1/3]

            q = -2* 1/(t[ip]-t[i]) *
            @SVector [-(I_tg[ip] - I_tg[i]) + t[ip]*(I_g[ip] - I_g[i]),
            (I_tg[ip] - I_tg[i]) - t[i]*(I_g[ip] - I_g[i])]

            r =  I_g2[ip] - I_g2[i]

            ℓ[i,ip] = QuadraticForm(P, q, r)
        end
    end

    return ℓ
end




"""
    ℓ = compute_transition_costs(g::AbstractArray)
Computes the transition costs `ℓ` given a discrete
function `g`.
"""
function compute_discrete_transition_costs(g::AbstractArray)
    T = Float64

    N = length(g)

    # Find sums of g, k*g, and g^2
    G1 = zeros(T, N)
    G2 = zeros(T, N)
    G3 = zeros(T, N)

    # The sums corresponding to transitioning from i to ip
    # i.e. not including the cost at ip
    for k=2:N
        G1[k] = G1[k-1] + g[k-1]
        G2[k] = G2[k-1] + (k-1)*g[k-1]
        G3[k] = G3[k-1] + g[k-1]^2
    end

    # The P-matrices only depend on the distance d=ip-i
    P_mats  = Vector{SMatrix{2,2,Float64,4}}(N-1)
    P_mats[1] = @SMatrix [1.0 0; 0 0]
    for d=2:N-1
        off_diag_elems = sum([k*(d - k) for k=0:d-1])
        P_mats[d] = @SMatrix [P_mats[d-1][1,1] + d^2    off_diag_elems;
        off_diag_elems            P_mats[d-1][1,1]]
    end

    P_mats = P_mats ./ (1.0:N-1).^2 # FIXME: Why can't this be done above in the loop?

    #P_invs = inv.(P_mats)

    ℓ = Array{QuadraticForm{T}}(N-1,N)

    for i=1:N-1
        for ip=i+1:N

            P = P_mats[ip-i]

            q = -2* 1/(ip-i) *
            @SVector [-(G2[ip] - G2[i]) + ip*(G1[ip] - G1[i]),
            (G2[ip] - G2[i]) - i*(G1[ip] - G1[i])]

            r =  G3[ip] - G3[i]

            ℓ[i,ip] = QuadraticForm(P, q, r)
        end
    end

    return ℓ
end
