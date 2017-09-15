using Base.Test

include(joinpath(Pkg.dir("DynamicApproximations"),"src","dev.jl"))

N = 50
t = linspace(0,4π,N)
g = sin

@time ℓ = dev.compute_transition_costs(g, t);
cost_last = dev.QuadraticPolynomial(0.0, 0.0, 0.0)

@time Λ = dev.find_optimal_fit(ℓ, cost_last, 7);

I_sols = [      [[1, 50]],
                [[1, 34, 50], [1, 17, 50]],
                [[1, 21, 30, 50]],
                [[1, 8, 19, 30, 50], [1, 21, 32, 43, 50]],
                [[1, 8, 19, 32, 43, 50]]
         ]

f_costs = [5.328255648628215, 4.8783936257642315, 1.7402065022125042, 0.9196880458290417, 0.09174442455649423]

@testset "Optimal Fit m=$m" for m = 1:5
    @time I, _, f = dev.recover_solution(Λ[m, 1], ℓ, cost_last)
    # Comparisons to solution found by brute force optimization
    @test I ∈ I_sols[m]
    @test f ≈ f_costs[m]    atol = 1e-8

    ## Test of brute force optimization
    @time I_bf, _, f_bf = dev.brute_force_optimization(ℓ, cost_last, m);
    @test I_bf ∈  I_sols[m]
end


##


@testset "Regularize ζ=$ζ" for ζ in 0.1:0.1:2

    Λ_reg = dev.regularize(ℓ, cost_last, ζ)

    # recover_optimal_index_set returns the cost inclusive the regularization penality,
    # revober optimal solution does not do so. It is arguably more interesting
    # to test cost including regularization.
    I, _, f_reg = dev.recover_optimal_index_set(Λ_reg[1])

    # Use the costs above to find out how many segments the solution should contain
    m_expected = indmin([f_costs[m] + m*ζ for m=1:5])

    @test m_expected == length(I) - 1
    @test I ∈ I_sols[m_expected]
    @test f_reg ≈ f_costs[m_expected] + ζ*m_expected
end
