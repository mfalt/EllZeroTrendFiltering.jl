using Base.Test, DynamicApproximations
using Interpolations

srand(31415)
## Test 1
# Compare the approximation costs computed using ℓ which is obtained from
# compute_discrete_transition_costs and computed via simple linear interpolation

function compute_cost(ℓ, I, Y)
    cost = 0.0;
    for k=1:length(I)-1
     cost += ℓ[I[k], I[k+1]](Y[k], Y[k+1])
    end
    return cost
end



t = linspace(0,π,50)
g = sin.(t)


ℓ = compute_discrete_transition_costs(g)

## Test 1a
I1 = [1,10,30,40,50]
Y1 = g[I1]

y1 = interpolate((I1,), Y1, Gridded(Linear()))[1:length(g)]

@test compute_cost(ℓ,I1,Y1) ≈ sum((y1[1:end-1]-g[1:end-1]).^2)

# Test with t grid
t_subsampled = [1, 10, 21, 33, 50]
ℓ_subsampled = compute_discrete_transition_costs(g, t=t_subsampled)

for i in 1:length(t_subsampled)-1
    for ip in (i+1):length(t_subsampled)
        t_i = t_subsampled[i]
        t_ip = t_subsampled[ip]
        @test ℓ_subsampled[i, ip] == ℓ[t_i, t_ip]
    end
end

## Test 1b
I2 = [1,15,30,36,50]
Y2 = randn(5)

y2 = interpolate((t[I2],), Y2, Gridded(Linear()))[t]

@test compute_cost(ℓ,I2,Y2) ≈ sum((y2[1:end-1]-g[1:end-1]).^2)

#plot(t, sin.(t))
#plot!(t[I1], Y1)
#plot!(t, y1, color="red")


## Test 2
# Check the transition costs for a simple example
ℓ = compute_discrete_transition_costs([1.0, 2.0, 3.0])

# (y - 1)^2
@test ℓ[1,2].P == [1.0 0.0; 0.0 0.0]
@test ℓ[1,2].q == [-2.0, 0]
@test ℓ[1,2].r == 1

# (y - 1)^2 + (y/2 + y'/2 - 2)^2
@test ℓ[1,3].P == [1.25 0.25; 0.25 0.25]
@test ℓ[1,3].q == [-4.0, -2.0]
@test ℓ[1,3].r == 5.0
