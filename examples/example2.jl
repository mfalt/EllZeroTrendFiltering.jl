using IterTools
using Plots

# The problem Σerror^2 + card(I)
# Heuristic

include(joinpath(Pkg.dir("DynamicApproximations"),"src","dev.jl"))
##

data = readdlm(joinpath(Pkg.dir("DynamicApproximations"),"examples","data","snp500.txt"))
data = data[1:300]

N = length(data)



@time ℓ = dev.compute_discrete_transition_costs(data);

#@time I, Y, f = dev.brute_force_optimization(ℓ, K-1);

cost_last = dev.QuadraticPolynomial(1.0, -2*data[end], data[end]^2)

start_time = time()
gc()
@time Λ = dev.find_optimal_fit(ℓ, cost_last, 10, 1.65);
println("Time: ", time()-start_time)

#println(dev.counter1, " ", dev.counter2)

for k=3:10
    println(k, " : ", dev.recover_optimal_index_set(Λ[k, 1], 1, N)[3])
end
###

@time I2, y2, f2 = dev.recover_optimal_index_set(Λ[7, 1], 1, N)
println(I2)
Y2, f2_2 = dev.find_optimal_y_values(ℓ, cost_last, I2)


I2, Y2, f3 = dev.recover_solution(Λ[7, 1], ℓ, cost_last, )
#println("Comparison: ", sqrt(f), " --- ", sqrt(f2))


#find_optimal_y_values

plot(layout=(1,2))
##
l = zeros(size(Λ))
for i=1:size(l,1)
    for j=1:size(l,2)
        if isassigned(Λ,i,j)
            l[i,j] = length(Λ[i, j])
        end
    end
end
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
