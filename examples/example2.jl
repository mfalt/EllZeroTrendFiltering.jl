using IterTools
using Plots
using EllZeroTrendFiltering
using DelimitedFiles

data = readdlm(joinpath(joinpath(dirname(@__FILE__),"data","snp500.txt")))
data = data[1:300]

N = length(data)

@time ℓ = compute_discrete_transition_costs(data);
cost_last = QuadraticPolynomial(1.0, -2*data[end], data[end]^2)

start_time = time()

@time Λ = pwq_dp_constrained(ℓ, cost_last, 10, 1.65)
println("Time: ", time()-start_time)

for k=3:10
    println(k, " : ", recover_optimal_index_set(Λ[k, 1])[3])
end
###



@time I2, y2, f2 = recover_optimal_index_set(Λ[7, 1], 1, N)
println(I2)
Y2, f2_2 = find_optimal_y_values(ℓ, cost_last, I2)


I2, Y2, f3 = recover_solution(Λ[7, 1], ℓ, cost_last, )
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

I2, y2, f2 = recover_solution(Λ[6, 1], 1, N)
println(I2)
Y2, _ = find_optimal_y_values(ℓ, I2)

plot(data)
plot!(I2, Y2)
