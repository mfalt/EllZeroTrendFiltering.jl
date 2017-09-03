using IterTools


include(joinpath(Pkg.dir("DynamicApproximations"),"src","dev.jl"))


data = readdlm(joinpath(Pkg.dir("DynamicApproximations"),"examples","data","snp500.txt"))
data = data[1:800]

N = length(data)



@time ℓ = dev.compute_discrete_transition_costs(data);



K = 4
#@time I, Y, f = dev.brute_force_optimization(ℓ, K-1);
###

cost_last = dev.QuadraticPolynomial(1.0, -2*data[end], data[end]^2)
Λ_0 = [dev.create_new_pwq(dev.minimize_wrt_x2_fast(ℓ[i, N], cost_last)) for i in 1:N-1];

start_time = time()
@profiler Λ = dev.find_optimal_fit(Λ_0, ℓ, 7);
println("Time: ", time()-start_time)
###

@time I2, y2, f2 = dev.recover_solution(Λ[7, 1], 1, N)
println(I2)
Y2, _ = dev.find_optimal_y_values(ℓ, I2)

#println("Comparison: ", sqrt(f), " --- ", sqrt(f2))


#find_optimal_y_values


l = zeros(size(Λ))
for i=1:size(l,1)
    for j=1:size(l,2)
        if isassigned(Λ,i,j)
            l[i,j] = length(Λ[i, j])
        end
    end
end



I2, y2, f2 = dev.recover_solution(Λ[6, 1], 1, N)
println(I2)
Y2, _ = dev.find_optimal_y_values(ℓ, I2)

plot(data)
plot!(I2, Y2)


plot(l')
