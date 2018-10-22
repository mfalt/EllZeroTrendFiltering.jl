using Test, EllZeroTrendFiltering
using Interpolations

Random.seed!(31415)

# Auxilliary function for evaluating the approximation error for a specific
# piecewise linear approximation (I, Y), given the transition costs ℓ
function compute_cost(ℓ, I, Y)
    cost = 0.0;
    for k=1:length(I)-1
     cost += ℓ[I[k], I[k+1]](Y[k], Y[k+1])
    end
    return cost
end


## Test 1
# Compare the approximation error for specific piecewise linear approximations
# (I, Y) when the costs are computed  by
# (1) using ℓ computed by compute_discrete_transition_costs and
#     then evaluated using the above auxilliary function
# (2) simple linear interpolation

g1 = sin.(range(0, stop=π, length=50))
I1 = [1,10,30,40,50]
Y1 = g1[I1]
t_ss1 = [1, 10, 21, 33, 50]

g2 = g1
I2 = [1,15,30,36,50]
Y2 = randn(5)
t_ss2 = t_ss1

g3 = randn(20)
I3 = [1; 3:2:19; 20]
Y3 = randn(length(I3))
t_ss3 = [1,5,10,15,20]

g4 = g3
I4 = 1:length(g3)
Y4 = g4
t_ss4 = t_ss3


for (g, I, Y, t_subsampled) in ((g1, I1, Y1, t_ss1), (g2, I2, Y2, t_ss2), (g3, I3, Y3, t_ss3))#, (g4, I4, Y4))
    ℓ = compute_discrete_transition_costs(g)
    cost1 = compute_cost(ℓ,I,Y)

    y = interpolate((I,), Y, Gridded(Linear()))(1:length(g))
    cost2 = sum((y[1:end-1]-g[1:end-1]).^2) # Note: cost at i=N should not be included

    @test cost1 ≈ cost2

    # Test transition cost with t grid
    ℓ_subsampled = compute_discrete_transition_costs(g, t_subsampled)

    for i in 1:length(t_subsampled)-1
        for ip in (i+1):length(t_subsampled)
            t_i = t_subsampled[i]
            t_ip = t_subsampled[ip]
            @test ℓ_subsampled[i, ip] == ℓ[t_i, t_ip]
        end
    end
end


## Test 2
# Check the transition costs for a very simple example
ℓ = compute_discrete_transition_costs([1.0, 2.0, 3.0])

# (y - 1)^2
@test ℓ[1,2].P == [1.0 0.0; 0.0 0.0]
@test ℓ[1,2].q == [-2.0, 0]
@test ℓ[1,2].r == 1

# (y - 1)^2 + (y/2 + y'/2 - 2)^2
@test ℓ[1,3].P == [1.25 0.25; 0.25 0.25]
@test ℓ[1,3].q == [-4.0, -2.0]
@test ℓ[1,3].r == 5.0
