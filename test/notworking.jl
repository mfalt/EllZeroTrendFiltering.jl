@everywhere include(joinpath(Pkg.dir("DynamicApproximations"),"src","dev.jl"))

@everywhere using Polynomials, IterTools, Plots

srand(9)
N = 20
ran = cumsum(0.1*randn(10*N))

t = linspace(0,1,N)

g6(t) = ran[floor(Int64, t*10(N-1)+1)]
ℓ = dev.compute_transition_costs(g6, t);

Λ_0 = [dev.create_new_pwq(dev.minimize_wrt_x2(ℓ[i, N], dev.QuadraticPolynomial{Float64}(0.,0.,0.))) for i in 1:N-1];
Λ, t2, _, _, _ = @timed dev.find_optimal_fit(Λ_0, ℓ, 18);
