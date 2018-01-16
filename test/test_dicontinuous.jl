function dp_discontinuous(C::Matrix{T}, M::Integer)
    N = size(C,2)
    D = Array{PiecewiseQuadratic{T}}(N, M)

    for m=2:M
        for i=1:N-m
            for ip=i+1:N-m+1
                D[i, m] = max(D[i, m], D[ip, m-1] + C[i, ip])
            end
        end
    end
end

##
compute_discrete_transition_costs(randn(srand(1), 30))
