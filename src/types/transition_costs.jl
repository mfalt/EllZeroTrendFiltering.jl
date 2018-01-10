## Transition cost object for continuous time approximation

struct TransitionCostContinuous{T}
    t::Vector{T}
    I_g::Vector{T}
    I_g2::Vector{T}
    I_tg::Vector{T}
end

function TransitionCostContinuous{T}(g, t::AbstractVector, tol=1e-3) where {T}
    N = length(t)
    I_g =  zeros(N)
    I_g2 = zeros(N)
    I_tg = zeros(N)
    for i=2:N
        I_g[i] = I_g[i-1] + quadgk(g, t[i-1], t[i], reltol=tol)[1]
        I_g2[i] = I_g2[i-1] + quadgk(τ -> g(τ)^2, t[i-1], t[i], reltol=tol)[1]
        I_tg[i] = I_tg[i-1] + quadgk(τ -> τ*g(τ), t[i-1], t[i], reltol=tol)[1]
    end
    TransitionCostContinuous{T}(t, I_g, I_g2, I_tg)
end


# getindex returns a quadratic form H([x1, x2]) that represents the
# transition cost from i to ip
function getindex(tc::TransitionCostContinuous, i::Integer, ip::Integer)
    t, I_g, I_g2, I_tg = tc.t, tc.I_g, tc.I_g2, tc.I_tg

    P = @SMatrix [(t[ip]-t[i])/3 (t[ip]-t[i])/6; (t[ip]-t[i])/6 (t[ip]-t[i])/3]

    q = -2* 1/(t[ip]-t[i]) *
    @SVector [-(I_tg[ip] - I_tg[i]) + t[ip]*(I_g[ip] - I_g[i]),
    (I_tg[ip] - I_tg[i]) - t[i]*(I_g[ip] - I_g[i])]

    r =  I_g2[ip] - I_g2[i]

    return QuadraticForm(P, q, r)
end

Base.size(tc::TransitionCostContinuous) = (length(tc.t)-1, length(tc.t))
Base.size(tc::TransitionCostContinuous, i::Integer) = size(tc)[i]

## Transition cost object for continuous time approximation

# Test non-uniform sampling
# non-uniform sampling needs to be handled

struct TransitionCostDiscrete{T}
    t::Vector{UInt32} # Should perhaps allow more general datatype
    P_matrices::Vector{SMatrix{2,2,T,4}}
    G1::Vector{T}
    G2::Vector{T}
    G3::Vector{T}
end

function TransitionCostDiscrete{T}(g::AbstractArray{T}, t=1:length(g)) where {T}
    N = length(g)

    # Find sums of g, k*g, and g^2
    G1 = zeros(T, N)
    G2 = zeros(T, N)
    G3 = zeros(T, N)
    # The sums corresponding to transitioning from [1, ip)
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

    P_mats = P_mats ./ (1.0:N-1).^2 # FIXME: Det var något problem här...
    TransitionCostDiscrete{T}(t, P_mats, G1, G2, G3)
end

# getindex returns a quadratic form H([x1, x2]) that represents the
# transition cost from i to ip
function getindex(tc::TransitionCostDiscrete, i::Integer, ip::Integer)

    q = -2* 1/(ip-i) *
    @SVector [-(tc.G2[ip] - tc.G2[i]) + ip*(tc.G1[ip] - tc.G1[i]),
    (tc.G2[ip] - tc.G2[i]) - i*(tc.G1[ip] - tc.G1[i])]

    r =  tc.G3[ip] - tc.G3[i]

    return QuadraticForm(tc.P_matrices[ip-i], q, r)
end

Base.size(tc::TransitionCostDiscrete) = (length(tc.t)-1, length(tc.t))
Base.size(tc::TransitionCostDiscrete, i::Integer) = size(tc)[i]
