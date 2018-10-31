using EllZeroTrendFiltering

N = 10
t = [0.0:N-1...]
I_g = rand(N)
I_g2 = rand(N)
I_tg = rand(N)

l = EllZeroTrendFiltering.TransitionCostContinuous{Float64}(t, I_g, I_g2, I_tg)

#l = compute_transition_costs(x -> exp(2^x), range(0, stop=2, length=N))

l = compute_transition_costs(x -> exp(x^2), range(0, stop=1, length=20))

V_N = QuadraticPolynomial(0.0, 0.0, 0.0)

M = 10

@time Λ =  pwq_dp_constrained(l, V_N, M)


unique_polys(Λ) = length(unique([λ.p for λ in Λ]))

sizes = [(isassigned(Λ,i,j) ? unique_polys(Λ[i,j]) : 0) for i = 1:size(Λ,1), j = 1:size(Λ,2)]
sizes = [(isassigned(Λ,i,j) ? length(Λ[i,j]) : 0) for i = 1:size(Λ,1), j = 1:size(Λ,2)]

# Check compared to bound R <= N
nom_sizes = fill(1, M, 1)*(20-1:-1:0)'

println(maximum(sizes - nom_sizes))
