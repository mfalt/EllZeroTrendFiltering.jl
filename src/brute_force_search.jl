"""
Brute force search for the optimal solution over all
index sets with m segemets. The costs are evaluated using least squares.

(Intended for verifying the dynamic programming algorithms)
"""
function brute_force_search(ℓ, V_0N::QuadraticPolynomial, m::Integer)
    cost_best = Inf

    I_best = []
    Y_best = []

    N = size(ℓ, 2)

    for I=IterTools.subsets(2:N-1, m-1)

        I = [1; I; N]
        P = zeros(m+1, m+1)
        q = zeros(m+1)
        r = 0

        # Add cost at right endpoint
        P[end,end] = V_0N.a
        q[end]     = V_0N.b
        r          = V_0N.c

        # Form quadratic cost function Y'*P*Y + q'*Y + r
        # corresponding to the y-values in the vector Y
        for j=1:m
            P[j:j+1,j:j+1] .+= ℓ[I[j], I[j+1]].P
            q[j:j+1] .+= ℓ[I[j], I[j+1]].q
            r += ℓ[I[j], I[j+1]].r
        end

        # find the optimal Y-vector, and compute the correspinding error
        Yopt = -(P \ q) / 2
        cost = Yopt' * P * Yopt + q' * Yopt + r

        if cost < cost_best
            cost_best = cost
            Y_best = Yopt
            I_best = I
        end
    end
    return I_best, Y_best, cost_best
end
