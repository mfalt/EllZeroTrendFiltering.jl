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
