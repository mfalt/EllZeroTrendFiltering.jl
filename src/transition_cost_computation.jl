
"""
    l = compute_transition_costs(g, t::AbstractArray)
Computes the transition costs `l` given a
function `g` and a time sequence `t`
"""
function compute_transition_costs(g, t::AbstractVector, tol=1e-3)
    T = Float64
    l_lazy = TransitionCostContinuous{T}(g, t)
    N = length(t)
    l = Array{QuadraticForm{T}}(undef, N-1,N)
    for i = 1:N-1
        for ip = i+1:N
            l[i,ip] = l_lazy[i,ip]
        end
    end
    # TODO move/add tests to Lazy, maybe only i to i+1
    for i = 1:N-1, ip = i+1:N
        P, q, r = l[i,ip].P, l[i,ip].q, l[i,ip].r
        minval = -q⋅(P\q)/4+r
        #TODO add epsilon on RHS?
        minval < 0 && error("Transition cost is negative for some values in `compute_transition_costs`, try decreasing rtol.")
    end

    return l
end




"""
    l = compute_transition_costs(g::AbstractArray; t=1:length(g))
Computes the transition costs `l` given a discrete
function `g`.
"""
function compute_discrete_transition_costs(g::AbstractArray, t=1:length(g))
    T = promote_type(eltype(g), Float64)
    l_lazy = TransitionCostDiscrete{T}(g, t)
    N = length(t)
    l = Array{QuadraticForm{T}}(undef, N-1,N)
    for i = 1:N-1
        for ip = i+1:N
            l[i,ip] = l_lazy[i,ip]
        end
    end
    return l
end
