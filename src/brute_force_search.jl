"""
Brute force search for the optimal solution over all
index sets with m segemets. The costs are evaluated using least squares.

(Intended for verifying the dynamic programming algorithms)
"""
function brute_force_search(l::AbstractTransitionCost{T}, V_N::QuadraticPolynomial{T}, m::Integer) where {T}
    cost_best = Inf
    I_best = Vector{Int64}(undef, m+1)
    Y_best = Vector{T}(undef, m+1)

    N = size(l, 2)

    I = zeros(Int64, m+1)
    I[1] = 1
    I[end] = N

    for I_inner=IterTools.subsets(2:N-1, m-1)

        I[2:m] .= I_inner

        P = SymTridiagonal(zeros(T, m+1), zeros(T, m+1))
        q = zeros(T, m+1)
        r = 0

        # Add cost at right endpoint
        P[end,end] = V_N.a
        q[end]     = V_N.b
        r          = V_N.c

        # Form quadratic cost function Y'*P*Y + q'*Y + r
        # corresponding to the y-values in the vector Y
        for j=1:m
            P.dv[j] += l[I[j], I[j+1]].P[1,1]
            P.dv[j+1] += l[I[j], I[j+1]].P[2,2]
            P.ev[j] += l[I[j], I[j+1]].P[1,2]
            q[j:j+1] .+= l[I[j], I[j+1]].q
            r += l[I[j], I[j+1]].r
        end

        # find the optimal Y-vector, and compute the correspinding error
        Yopt = -(P \ q) / 2
        cost = dot(Yopt, P*Yopt) + q' * Yopt + r

        if cost < cost_best
            cost_best = cost
            Y_best .= Yopt
            I_best .= I
        end
    end
    return I_best, Y_best, cost_best
end




"""
Brute force search for best subset selection.
I.e. minimize    ||Ax - b||_2
     subject to  card(x[c:end]) <= m
"""
function brute_force_search(A::Matrix{T}, b::Vector{T}, m::Integer, c::Integer=0) where {T}
    cost_best = Inf
    I_best = Vector{Int64}(undef, m+c)
    u_best = Vector{T}(undef, m+c)

    N = size(A, 2)

    I = zeros(Int64, m+c)

    for I_inner=IterTools.subsets(c+1:N, m)

        I .= [1:c; I_inner]

        # find the optimal u_opt vector
        u_opt = A[:, I] \ b
        cost = norm(A[:, I]*u_opt - b)^2 # Squared 2-norm is used everywhere else

        #println("$I_inner : $cost")

        if cost < cost_best
            cost_best = cost
            u_best .= u_opt
            I_best .= I
        end
    end
    return I_best, u_best, cost_best
end



# Either initial conditions are free or they are zero
function generate_markov_matrix(sys::ControlSystems.StateSpace, N; free_intial_conditions=false)
    y, _ = ControlSystems.impulse(sys, N-1)
    T = MatrixDepot.matrixdepot("toeplitz", y[1:N],  zeros(N))

    if free_intial_conditions
        sys_initial = ControlSystems.ss(sys.A, sys.B, Matrix(I, 2, 2), 0, 1)
        y_initial, _ = ControlSystems.impulse(sys_initial, N)
        return [y_initial[2:end, :] T]
    else
        return T
    end
end
