# PiecewiseLinearContinuousFitting
# Piecewise Linear Continouous Interpolation Constraints
# Sparse L2 Optimal Fitting subject Continuity Constraints
# Integrated Square

include("datatypes.jl")

using Polynomials
using IterTools

t = linspace(0,1,30)




g = Poly([0,0,1,-1])

ℓ = dev.compute_transition_costs(g, t)



[I, Y] = brute_force_optimization(ℓ, K)
  

x = linspace(0,1,100)

#closefig()

using PyPlot
figure(1)
plot(x, g(x), "g", t[I], Y, "r")
