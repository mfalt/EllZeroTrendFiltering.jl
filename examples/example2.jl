using IterTools
using Plots

include(joinpath(Pkg.dir("DynamicApproximations"),"src","dev.jl"))
##

data = readdlm(joinpath(Pkg.dir("DynamicApproximations"),"examples","data","snp500.txt"))
data = data[1:300]

N = length(data)



@time ℓ = dev.compute_discrete_transition_costs(data);

#@time I, Y, f = dev.brute_force_optimization(ℓ, K-1);

cost_last = dev.QuadraticPolynomial(1.0, -2*data[end], data[end]^2)
Λ_0 = [dev.create_new_pwq(dev.minimize_wrt_x2(ℓ[i, N], cost_last)) for i in 1:N-1];
start_time = time()
gc()
@time Λ = dev.find_optimal_fit(Λ_0, ℓ, 10, Inf);
println("Time: ", time()-start_time)

#println(dev.counter1, " ", dev.counter2)

for k=3:10
    println(k, " : ", dev.recover_solution(Λ[k, 1], 1, N)[3])
end
###

@time I2, y2, f2 = dev.recover_solution(Λ[7, 1], 1, N)
println(I2)
Y2, _ = dev.find_optimal_y_values(ℓ, I2)

#println("Comparison: ", sqrt(f), " --- ", sqrt(f2))


#find_optimal_y_values

plot(layout=(1,2))
##
l = [(isassigned(Λ,i,j) ? length(Λ[i,j]) : 0) for i = 1:size(Λ,1), j = 1:size(Λ,2)]
sum(l)
##

plot!(l', subplot=2)

sum(l)
##

I2, y2, f2 = dev.recover_solution(Λ[6, 1], 1, N)
println(I2)
Y2, _ = dev.find_optimal_y_values(ℓ, I2)

plot(data)
plot!(I2, Y2)
