using Base.Test

include(joinpath(Pkg.dir("DynamicApproximations"),"src","dev.jl"))

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

# (3x+1)^2 + 1 = 9x^2 + 6x + 2
@test dev.find_minimum(dev.QuadraticPolynomial(9, 6, 2)) == (-1/3, 1)

# (3x+1)^2 + 1 = 9x^2 + 6x + 2
@test dev.unsafe_minimum(dev.QuadraticPolynomial(9.0, 6, 2)) == 1
