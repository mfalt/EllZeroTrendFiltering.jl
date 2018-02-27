using IterTools
using Plots

# The problem Σerror^2 + card(I)
# Heuristic

include(joinpath(Pkg.dir("EllZeroTrendFiltering"),"src","EllZeroTrendFiltering.jl"))
##
DA = EllZeroTrendFiltering

data = readdlm(joinpath(Pkg.dir("EllZeroTrendFiltering"),"examples","data","snp500.txt"))
data = data[1:300]

N = length(data)



@time ℓ = DA.compute_discrete_transition_costs(data);

#@time I, Y, f = brute_force_search(ℓ, K-1);

cost_last = DA.QuadraticPolynomial(1.0, -2*data[end], data[end]^2)

start_time = time()
gc()
@time Λ = DA.pwq_dp_constrained(ℓ, cost_last, 10, 1.65);
println("Time: ", time()-start_time)

#println(counter1, " ", counter2)

for k=3:10
    println(k, " : ", DA.recover_optimal_index_set(Λ[k, 1])[3])
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
