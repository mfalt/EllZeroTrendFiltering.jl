using Base.Test

include("../src/dev.jl")

# Subtraction of quadratic polynomials
@test dev.QuadraticPolynomial(4.0, -2.0, 3.0) - dev.QuadraticPolynomial(1.0, 1.0, 2.0) == dev.QuadraticPolynomial(3.0,-3.0,1.0)


### Minimization of quadratic polynomials ###

@test dev.find_minimum(dev.QuadraticPolynomial(1 ,0 ,0)) == (0, 0)
@test dev.find_minimum(dev.QuadraticPolynomial(1 ,0 ,1)) == (0, 1)

# 2*x*(x+1) + 1
@test dev.find_minimum(dev.QuadraticPolynomial(2, 2, 1)) == (-0.5, 0.5)

# (x+1)^2 + 1
@test dev.find_minimum(dev.QuadraticPolynomial(1, 2, 2)) == (-1, 1)

# 2*((x+1)^2 + 1)
@test dev.find_minimum(dev.QuadraticPolynomial(2, 4, 4)) == (-1, 2)
