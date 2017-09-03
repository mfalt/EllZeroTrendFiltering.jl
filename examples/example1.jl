using Plots
using Convex
using Mosek
using SCS



data = readdlm(joinpath(Pkg.dir("DynamicApproximations"),"examples","data","snp500.txt"))
n = size(data)[1]

#p = plot(1:length(data), data)

N = length(data)
H =  spdiagm((ones(N-2), -2*ones(N-2), ones(N-2)), (0,1,2))


x = Variable(N)

λ = 50
problem = minimize(0.5*sumsquares(x - data) + λ*norm(H*x, 1))

# Solve the problem by calling solve!
solve!(problem, MosekSolver())

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimal value
problem.optval

plot(data)
plot!(evaluate(x))
