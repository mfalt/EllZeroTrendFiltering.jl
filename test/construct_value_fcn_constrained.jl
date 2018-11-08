using Test, EllZeroTrendFiltering

t = range(0, stop=4π, length=50)
g = sin

@time l, χ, V_N = compute_problem_data_pwl(g, t; tol=1e-10)
V_N = QuadraticPolynomial(0.0, 0.0, 0.0)

@time Λ =  construct_value_fcn_constrained(l, χ, V_N, 7);

I_sols = [      [[1, 50]],
                [[1, 34, 50], [1, 17, 50]],
                [[1, 21, 30, 50]],
                [[1, 8, 19, 30, 50], [1, 21, 32, 43, 50]],
                [[1, 8, 19, 32, 43, 50]]
         ]

f_costs = [5.328255648628215, 4.8783936257642315, 1.7402065022125042, 0.9196880458290417, 0.09174442455649423]


@testset "Optimal Fit m=$m" for m = 1:5
    I, _, f = recover_solution(Λ[m, 1], l, V_N)

    @test I ∈ I_sols[m]
    @test f ≈ f_costs[m]    atol = 1e-8
end

##
@testset "Brute force m=$m" for m = 1:4
    I_bf, _, f_bf = brute_force_search(l, V_N, m);

    @test I_bf ∈  I_sols[m]
    @test f_bf ≈ f_costs[m]    atol = 1e-8
end

##


@testset "Regularize ζ=$ζ" for ζ in 0.1:0.1:2

    Λ_reg = construct_value_fcn_regularized(l, χ, V_N, ζ)

    # recover_optimal_index_set returns the cost inclusive the regularization penality,
    # revober optimal solution does not do so. It is arguably more interesting
    # to test cost including regularization.
    I, _, f_reg = recover_optimal_index_set(Λ_reg[1])

    # Use the costs above to find out how many segments the solution should contain
    m_expected = argmin([f_costs[m] + m*ζ for m=1:5])

    @test m_expected == length(I) - 1
    @test I ∈ I_sols[m_expected]
    @test f_reg ≈ f_costs[m_expected] + ζ*m_expected
end