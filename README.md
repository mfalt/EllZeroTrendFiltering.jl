# EllZeroTrendFiltering

[![Build Status](https://travis-ci.com/mfalt/EllZeroTrendFiltering.jl.svg?token=a1HpLsx1pmUnusz71XN8&branch=master)](https://travis-ci.com/mfalt/EllZeroTrendFiltering.jl)

[![codecov](https://codecov.io/gh/mfalt/EllZeroTrendFiltering.jl/branch/master/graph/badge.svg?token=nt4j2gNB2l)](https://codecov.io/gh/mfalt/EllZeroTrendFiltering.jl)


This package solves the problem of piecewise linear, *continuous*, approximation subject to either a hard limit or a regularization penalty on the number of break points. An exact solution is obtained using dynamic program over piecewise quadratic function, which avoids the combinatorial complexity of a naive approach.

### Mathematical Description

We want to find a piecewise linear, continuous function ![f_{I,Y}](figures/f_IY.svg) with few segments, where `I` is the set of breakpoints, and `Y` are the values of the function at those breakpoints. This problem can be formulated as a constrained problem

![Problem Formulation Constrained](figures/problemConstrained.svg)

where `M` is the number of segments, or as a regulartized problem

![Problem Formulation Regularized](figures/problemRegularized.svg)

These problems can be solved for a discrete set of points `g` using
```julia
I, Y, v = fit_pwl_constrained(g, M)
I, Y, v = fit_pwl_regularized(g, ζ)
```
where the resulting function satisfies `f(I[k]) = Y[k]`,
or for a general function `g: ℝ ⟶ ℝ` with
```julia
I, Y, v = fit_pwl_constrained(g, t, M)
I, Y, v = fit_pwl_regularized(g, t, ζ)
```
where `t` is the set of possible breakpoints and the resulting function satisfies `f(t[I[k]]) = Y[k]`.
### Example: Constrained approximation

Find best continouous piecewise approximations with up to M segments
```julia
using EllZeroTrendFiltering, Plots
#Get sample data
N = 400
data = snp500_data()[1:N]

#Find solutions to the constrained problem, with up to M=10 segments
M = 10
Ivec, Yvec, fvec = fit_pwl_constrained(data, M)

#Plot original data
plot(data, l=:black, lab = "SNP500")
#Plot solution with 5 segments
plot!(Ivec[5], Yvec[5], l=2, m=:circle, lab="m=5, cost = $(round(fvec[5],digits=3))")
#Plot solution with 10 segments
plot!(Ivec[M], Yvec[M], l=2, m=:circle, lab="m=$M, cost = $(round(fvec[M],digits=3))")

```
![Example figure](figures/snp500.svg)

### Example: Regularization

Find best continouous piecewise approximations with cost ζ per breakpoint.
```julia
using EllZeroTrendFiltering, Plots

g(x) = sin(x) + 0.5sin(3.5x) + 0.5sin(5.1x)
t = range(0, stop=2π, length=201)

plot(g, t, l=(2,:black), lab="sin(x) + 0.5sin(3.5x) + 0.5sin(5.1x)")
for ζ ∈ [0.1, 0.002]
    # Minimize ∫(f(x)-g(x))²dx + ζ⋅||d²f/dx²||₀
    # Will automatically integrate the function to compute the costs
    I, Y, cost = fit_pwl_regularized(g, t, ζ)

    plot!(t[I], Y, l=2, m=:circle, lab = "l2-norm=$(round(cost,digits=3)), zeta=$ζ")
end
plot!() # Show plot
```
![Example figure](figures/sin.svg)
