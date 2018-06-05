using Plots
using Convex
using Mosek
using SCS
using Interpolations
using DelimitedFiles

# This is a julia implementaiton of an example
# from "ℓ1 Trend Filtering" by Kim et al (2009)

I_sets = [
[1, 2000],
[1, 968, 2000],
[1, 339, 897, 2000],
[1, 337, 983, 1145, 2000],
[1, 317, 791, 885, 1207, 2000],
[1, 350, 636, 755, 884, 1208, 2000],
[1, 349, 631, 788, 834, 991, 1150, 2000],
[1, 349, 631, 787, 835, 988, 1191, 1853, 2000],
[1, 350, 636, 755, 887, 922, 976, 1181, 1853, 2000],
[1, 365, 515, 545, 628, 769, 837, 987, 1190, 1853, 2000]
]

cost_l0 = [find_optimal_y_values(ℓ, I_sets[k])[2] for k=1:10]


#---

function is_local_max(y)
    v = zeros(y, Bool)
    for k=2:length(y)-1
        v[k] = (y[k] >= y[k-1]) && (y[k] >= y[k+1])
    end
    return v
end

#---
data = readdlm(joinpath(joinpath(dirname(@__FILE__),"data",,"snp500.txt"))

N = length(data)
H =  spdiagm((ones(N-2), -2*ones(N-2), ones(N-2)), (0,1,2))

x = Variable(N)


#---
r = 6
I = I_sets[r]
Y, _ = find_optimal_y_values(ℓ, I)
y_l0 = interpolate((I,), Y, Gridded(Linear()))[1.0:N]





λ = 100
problem = minimize(0.5*sumsquares(x - data) + λ*norm(H*x, 1))


solve!(problem, MosekSolver([("QUIET", true)]))
println("Opt. val: ",problem.optval, "   (problem status: ", problem.status, ")")

y_l1 = evaluate(x)[:]



v = sortperm(abs.(H*y_l1) + 1*is_local_max(abs.(H*y_l1)))

plot(H*y_l1)
plot!(v[end-r+2:end], zeros(r-1), m=:circle)

I_l1 = [1; sort(v[end-r+2:end]); N]
Y_l1, f_l1 = find_optimal_y_values(ℓ, I_l1)

y_l10 = interpolate((I_l1,), Y_l1, Gridded(Linear()))[1.0:N]

println("ℓ1: ", sum((evaluate(x) - data).^2))
println("ℓ0/ℓ1: ", sum((y_l10 - data).^2))
println("ℓ0: ", sum((y_l0 - data).^2))


#subplot(211)
plot(data)
#plot!(evaluate(x), color="green", linewidth=3)
plot!(I, Y, color="red", linewidth=3)
plot!(I_l1, Y_l1, l=(:cyan, :dash), linewidth=3)
#---

λ_vec = 10 .^ range(1, stop=4, length=60)
cost = zeros(length(λ_vec), 10)

for (k, λ)=enumerate(λ_vec)
    problem = minimize(0.5*sumsquares(x - data) + λ*norm(H*x, 1))

    solve!(problem, MosekSolver([("QUIET", true)]))
    println("Opt. val: ",problem.optval, "   (problem status: ", problem.status, ")")

    y_l1 = evaluate(x)[:]
    v = sortperm(abs.(H*y_l1) + 1*is_local_max(abs.(H*y_l1)))
    for r=1:10
        I_l1 = [1; sort(v[end-r+2:end]); N]
        _, f_l1 = find_optimal_y_values(ℓ, I_l1)
        cost[k, r] = f_l1
    end
end

#---

#---
plot(log10(λ_vec), cost[:,3:end], show=true, yaxis=((0,10), 0:0.5:10) )

plot!([1.5, 2.5], [1,1]*cost_l0[3:end]')

#---
[minimum(cost,1)' cost_l0]


#---
 costs = [42.8484   42.8484
  9.25796   9.25796
  4.46404   4.26137
  3.13439   3.02902
  3.09011   2.71251
  2.68868   2.27944
  2.16277   2.03848
  1.94739   1.8517
  1.89786   1.76981
  1.75176   1.59654]

plot(3:10,costs[3:end,:])
