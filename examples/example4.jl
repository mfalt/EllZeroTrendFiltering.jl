# A comparison between ℓ_0 and ℓ_1 trend filtering

using Plots
using Convex
using Mosek
using EllZeroTrendFiltering

g = snp500_data()[1:300]


N = length(g)


## ℓ_1 trend filtering, see "ℓ_1 Trend Filtering", Kim et al. (2009)
H =  spdiagm(0=>ones(N-2), 1=>-2*ones(N-2), 2=>ones(N-2))

x = Variable(N)

λ = 2.0
problem = minimize(sumsquares(x - g) + λ * norm(H*x, 1))

solve!(problem, Mosek.Optimizer)
#solve!(problem, SCS.Optimizer)
println(problem.status)

x_l1 = evaluate(x)[:]
f_l1 = sum((x_l1 - g).^2)


## ℓ_0 trend filtering
ζ = 0.02
I_l0, Y_l0, f_l0 = fit_pwl_regularized(g, ζ)



plot(g, label="SP500", size=(1000,550));
plot!(x_l1, lw=2.5, label="l1 reg. lambda = $λ");
plot!(I_l0, Y_l0, lw=2.5, label="l0 reg. zeta = $ζ")

# plot(H*x.value);
# plot!(I_l0[2:end], 0.1*diff(Y_l0) ./ diff(I_l0), m=:o, lw=0)
#
