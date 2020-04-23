using Plots
using Convex
using Mosek, MosekTools
using SCS
using DelimitedFiles
using SparseArrays

data = readdlm(joinpath(dirname(@__FILE__),"data","snp500.txt"))
n = size(data)[1]

#p = plot(1:length(data), data)

N = length(data)
H =  spdiagm(0=>ones(N-2), 1=>-2*ones(N-2), 2=>ones(N-2))


x = Variable(N)

λ = 50
problem = minimize(0.5*sumsquares(x - data) + λ*norm(H*x, 1))

# Solve the problem by calling solve!
solve!(problem, Mosek.Optimizer)
solve!(problem, SCS.Optimizer)

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimal value
problem.optval

plot(data)
plot!(evaluate(x))
