using Base.Test, IterTools
using DynamicApproximations: find_optimal_y_values

data = readdlm(joinpath(Pkg.dir("DynamicApproximations"),"examples","data","snp500.txt"))
N = 300
data = data[1:N]

@time ℓ = compute_discrete_transition_costs(data);


cost_last = QuadraticPolynomial(1.0, -2*data[N], data[N]^2)
Λ = find_optimal_fit(ℓ, cost_last, 10, Inf);

@testset "Discrete equivalent recover_optimal_index_set and find_optimal_y_values m=$m" for m in 1:10
    I2, y2, f2 = recover_optimal_index_set(Λ[m, 1])
    Y2, f2_2 = find_optimal_y_values(ℓ, cost_last, I2)

    @test Y2[1] ≈ y2  rtol=1e-10
    @test f2 ≈ f2_2   rtol=1e-10
end
#Test that recover_optimal_index_set gives same cost as find_optimal_y_values

@testset "Discrete compare to bruteforce m=$m" for m in 1:3
    I2, y2, f2 = recover_optimal_index_set(Λ[m, 1])
    Y2, f2_2 = find_optimal_y_values(ℓ, cost_last, I2)
    I3, Y3, f3 = brute_force_optimization(ℓ, cost_last, m)

    @test I2 == I3
    @test Y2 ≈ Y3    rtol=1e-10
    @test f2 ≈ f3    rtol=1e-10
    @test f2_2 ≈ f3  rtol=1e-10
end
