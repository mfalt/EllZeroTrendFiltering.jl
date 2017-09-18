
data = readdlm(joinpath(Pkg.dir("DynamicApproximations"),"examples","data","snp500.txt"))
N = 2000



cost_last = dev.QuadraticPolynomial(1.0, -2*data[N], data[N]^2)

Ns = vcat(collect.([2:100; 100:10:1000; 1100:100:2000])...)
times = Array{Float64}(length(Ns))

ζ = 1.0
for (i,N) in enumerate(Ns)
    println(i,",",N)
    data2 = data[1:N]
    gc()
    @time ℓ = dev.compute_discrete_transition_costs(data2);

    Λ_reg, t, _, _, _ = @timed dev.regularize(ℓ, cost_last, ζ)
    times[i] = t
end

plot(Ns, times)
plot!(Ns, Ns.^3/Ns[end]^3*times[end])
