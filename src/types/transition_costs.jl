## Transition cost object for continuous time approximation

struct TransitionCostContinuous{T}
    t::Vector{T}
    I_g::Vector{T}
    I_g2::Vector{T}
    I_tg::Vector{T}
end

function TransitionCostContinuous{T}(g, t::AbstractVector, tol=1e-3) where {T}
    N = length(t)
    I_g =  fill(zero(T), N)
    I_g2 = fill(zero(T), N)
    I_tg = fill(zero(T), N)
    for i=2:N
        I_g[i] = I_g[i-1] + quadgk(g, t[i-1], t[i], rtol=tol)[1]
        I_g2[i] = I_g2[i-1] + quadgk(τ -> g(τ)^2, t[i-1], t[i], rtol=tol)[1]
        I_tg[i] = I_tg[i-1] + quadgk(τ -> τ*g(τ), t[i-1], t[i], rtol=tol)[1]
    end
    l = TransitionCostContinuous{T}(t, I_g, I_g2, I_tg)
    #_test_positivity(l)
    return l
end

""" Check that all transition costs are positive for all x ∈ ℝ^d.
    This may not be the case due to numerical errors in their computation.
    If they are negative for some x, it causes the optimization to break down.
"""
function _test_positivity(l::TransitionCostContinuous)
    N = size(l, 2)
    for i=1:N-1, ip=i+1:min(i+2,N) # It should be sufficient to test only for elements on the main diagonal, one more diagonal is considered just in case
        #QUESTION: add epsilon on RHS?
        find_minimum(l[i,ip])[2] < 0 && error("Transition cost is negative for some values in `compute_transition_costs`, try decreasing rtol.")
    end
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

struct TransitionCostDiscrete{T,TimeType}
    t::TimeType
    P_matrices::Vector{SMatrix{2,2,T,4}}
    G1::Vector{T}
    G2::Vector{T}
    G3::Vector{T}
end

function TransitionCostDiscrete{T}(g::AbstractArray, t::TimeType=1:length(g)) where {T, TimeType<:AbstractArray{<:Integer}}
    N = length(g)

    if t[1] != 1 || t[end] != N
        @warn "TransitionCostDiscrete: The supplied grid t only covers the range ($(t[1]),$(t[end])) while the range of indices for g is (1,$(length(g))). No costs will be considered outside the range of t."
    end

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
    P_matrices  = Vector{SMatrix{2,2,T,4}}(undef, N-1)
    P_matrices[1] = @SMatrix [1.0 0; 0 0]
    for d=2:N-1
        off_diag_elems = sum([k*(d - k) for k=0:d-1])
        P_matrices[d] = @SMatrix [P_matrices[d-1][1,1] + d^2    off_diag_elems;
        off_diag_elems            P_matrices[d-1][1,1]]
    end

    P_matrices = P_matrices ./ (1.0:N-1).^2 # FIXME: Some problem here
    TransitionCostDiscrete{T,TimeType}(t, P_matrices, G1, G2, G3)
end

# getindex returns a quadratic form H([x1, x2]) that represents the
# transition cost from i to ip
function getindex(tc::TransitionCostDiscrete, i::Integer, ip::Integer)
    t_i = tc.t[i]
    t_ip = tc.t[ip]
    q = -2* 1/(t_ip-t_i) *
    @SVector [-(tc.G2[t_ip] - tc.G2[t_i]) + t_ip*(tc.G1[t_ip] - tc.G1[t_i]),
    (tc.G2[t_ip] - tc.G2[t_i]) - t_i*(tc.G1[t_ip] - tc.G1[t_i])]

    r =  tc.G3[t_ip] - tc.G3[t_i]

    return QuadraticForm(tc.P_matrices[t_ip-t_i], q, r)
end


Base.size(tc::TransitionCostDiscrete) = (length(tc.t)-1, length(tc.t))
Base.size(tc::TransitionCostDiscrete, i::Integer) = size(tc)[i]
