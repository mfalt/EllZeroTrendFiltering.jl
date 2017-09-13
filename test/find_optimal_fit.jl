using Base.Test

include(joinpath(Pkg.dir("DynamicApproximations"),"src","dev.jl"))

N = 50
t = linspace(0,4π,N)
g = sin

@time ℓ = dev.compute_transition_costs(g, t);
cost_last = dev.QuadraticPolynomial(0.0, 0.0, 0.0)
#Λ_0 = [dev.create_new_pwq(dev.minimize_wrt_x2(ℓ[i, N], ) for i in 1:N-1];

@time Λ = dev.find_optimal_fit(ℓ, cost_last, 7);

@time I_dp3, _, f_dp3 = dev.recover_solution(Λ[3, 1], ℓ, cost_last)
@time I_dp4, _, f_dp4 = dev.recover_solution(Λ[4, 1], ℓ, cost_last)
@time I_dp5, _, f_dp5 = dev.recover_solution(Λ[5, 1], ℓ, cost_last)

@testset "Optimal Fit" begin
    # Comparisons to solutions found by brute force optimization
    @test I_dp3 ==  [1, 21, 30, 50]
    @test I_dp4 ∈ ([1, 8, 19, 30, 50], [1, 21, 32, 43, 50]) # Two symmetric solutions
    @test I_dp5 == [1, 8, 19, 32, 43, 50]

    @test f_dp3 ≈ 1.7402065022125042
    @test f_dp4 ≈ 0.9196880458290417
    @test f_dp5 ≈ 0.09174442455649423

    ## Test of brute force optimization
    m = 3
    @time I_bf, _, f_bf = dev.brute_force_optimization(ℓ, m-1);
    @test I_bf ==  [1, 21, 30, 50]

end
