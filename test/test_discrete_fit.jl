using Base.Test
using DynamicApproximations: find_optimal_y_values


for problem_fcn in ["straight_line_problem", "square_wave_problem", "exp_problem", "snp500_problem"]

    include(joinpath(Pkg.dir("DynamicApproximations"),"test","problems", problem_fcn * ".jl"))
##
    g, ζvec, I_sols, f_sols = @eval $(Symbol(problem_fcn))()

    M = length(I_sols)


    @time ℓ = compute_discrete_transition_costs(g);
    cost_last = QuadraticPolynomial(1.0, -2*g[end], g[end]^2)

    Λ = pwq_dp_constrained(ℓ, cost_last, M, Inf);

    @testset "Data set: $problem_fcn, constrained with m=$m" for m in 1:M
        I, Y, f = recover_solution(Λ[m, 1], ℓ, cost_last)

        @test isempty(I_sols[m]) || I == I_sols[m]
        @test f ≈ f_sols[m]   rtol=1e-10 atol=1e-10
    end

    @testset "Data set: $problem_fcn, regularization with ζ=$ζ" for ζ in ζvec

        Λ_reg = pwq_dp_regularized(ℓ, cost_last, ζ)

        # the cost f returned from recover_optimal_index_set includes the
        # regularization penality mζ
        I, _, f_reg = recover_optimal_index_set(Λ_reg[1])

        # Use the costs above to find out how many segments the solution should contain
        m_expected = indmin([f_sols[m] + m*ζ for m=1:length(f_sols)])

        @test m_expected == length(I) - 1
        @test isempty(I_sols[m_expected]) || I == I_sols[m_expected]
        @test f_reg ≈ f_sols[m_expected] + ζ*m_expected
        
    end

end
