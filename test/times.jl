using Plots
plotly()

data = readdlm(joinpath(Pkg.dir("EllZeroTrendFiltering"),"examples","data","snp500.txt"))
N = 2000


Ns = vcat(collect.([2:100; 100:10:1000; 1100:100:2000])...)
times = Array{Float64}(length(Ns))

N = 2000
M = 15
data2 = data[1:N]
ζ = 1.
ℓ = compute_discrete_transition_costs(data2);
cost_last = QuadraticPolynomial(1.0, -2*data[N], data[N]^2)
@time Λ = pwq_dp_constrained(ℓ, cost_last, M, 1.3);
I, _, f = recover_solution(Λ[M, 1], ℓ, cost_last)
#@time Λ_reg = pwq_dp_regularized(ℓ, cost_last, ζ, 7.1)
I, _, f_reg = recover_optimal_index_set(Λ_reg[1])
find_minimum(Λ_reg[1])
lengths3 = [(isassigned(Λ, i, j) ? length(Λ[i,j]) : 0) for i = 1:size(Λ,1), j = 1:size(Λ,2)]
maximum(lengths3)
mean(lengths3)
ζ = 1.0
for (i,N) in enumerate(Ns)
    println(i,",",N)
    data2 = data[1:N]
    gc()
    @time ℓ = compute_discrete_transition_costs(data2);
    cost_last = QuadraticPolynomial(1.0, -2*data[N], data[N]^2)
    Λ_reg, t, _, _, _ = @timed pwq_dp_regularized(ℓ, cost_last, ζ, 7.03)
    times[i] = t
end
(i,N) = 1, 2000

plot(Ns, times)
plot!(Ns, Ns.^3/Ns[end]^3*times[end])


plotly()
x, y, x_all, y_all = get_vals(Λ[5,1])
plot(x,y,l=(2,:green), lab="$i,");
for (xi,yi) in zip(x_all,y_all);
    plot!(xi, yi, l=:black, lab="");
    plot!([xi[1]],[yi[1]], m=(:cross,:purple), lab="");
end
plot!()
remove_over(Λ_reg[601], 6.)

# # 9.65 - 9.75s on julia-903644385b with armv7l-unknown-linux-gnueabihf
# t1 = time()
# println("A: $(randn(10,10))")
# println(time()-t1)
