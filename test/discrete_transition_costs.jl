using Base.Test
using PyPlot
using Interpolations
include("../src/dev.jl")

t = linspace(0,π,50)
y = sin.(t)


ℓ = dev.compute_discrete_transition_costs(y)


I = [1,10,30,40,50]
Y = sin.(t[I])

cost = 0.0;
for k=1:length(I)-1
 i = I[k]
 ip = I[k+1]
 println(i, ":", ip)
 cost += ℓ[i, ip](Y[k], Y[k+1])
end

y2 = interpolate((t[I],), y[I], Gridded(Linear()))[t]

plot(t, sin.(t))
plot!(t[I], Y)
plot!(t, y2, "r--")
println(cost)

println("Error ℓ: ", cost, "    Error interp: ", sum((y-y2).^2))



##

ℓ = dev.compute_discrete_transition_costs([1.0, 2.0, 3.0])

# (y - 1)^2
@test ℓ[1,2].P == [1.0 0.0; 0.0 0.0]
@test ℓ[1,2].q == [-2.0, 0]
@test ℓ[1,2].r == 1

# (y - 1)^2 + (y/2 + y'/2 - 2)^2
@test ℓ[1,3].P == [1.25 0.25; 0.25 0.25]
@test ℓ[1,3].q == [-4.0, -2.0]
@test ℓ[1,3].r == 5.0
