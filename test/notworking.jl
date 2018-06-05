# TODO What is this?
#@everywhere include(joinpath(dirname(@__FILE__),"..","src","jl"))

@everywhere using Polynomials, IterTools, Plots

srand(9)
N = 20
ran = cumsum(0.1*randn(10*N))

t = range(0, stop=1, length=N)

g6(t) = ran[floor(Int64, t*10(N-1)+1)]
ℓ = compute_transition_costs(g6, t);

Λ_0 = [create_new_pwq(EllZeroTrendFiltering.minimize_wrt_x2(ℓ[i, N], QuadraticPolynomial{Float64}(0.,0.,0.))) for i in 1:N-1];
Λ, t2, _, _, _ = @timed pwq_dp_constrained(Λ_0, ℓ, 18);
