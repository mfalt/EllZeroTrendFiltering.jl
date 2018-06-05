using Test
using EllZeroTrendFiltering: unsafe_minimum

# Subtraction of quadratic polynomials
@test QuadraticPolynomial(4.0, -2.0, 3.0) - QuadraticPolynomial(1.0, 1.0, 2.0) == QuadraticPolynomial(3.0,-3.0,1.0)


### Minimization of quadratic polynomials ###
@test find_minimum(QuadraticPolynomial(1 ,0 ,0)) == (0, 0)
@test find_minimum(QuadraticPolynomial(1 ,0 ,1)) == (0, 1)

# 2*x*(x+1) + 1
@test find_minimum(QuadraticPolynomial(2, 2, 1)) == (-0.5, 0.5)

# (x+1)^2 + 1
@test find_minimum(QuadraticPolynomial(1, 2, 2)) == (-1, 1)

# 2*((x+1)^2 + 1)
@test find_minimum(QuadraticPolynomial(2, 4, 4)) == (-1, 2)

# (3x+1)^2 + 1 = 9x^2 + 6x + 2
@test find_minimum(QuadraticPolynomial(9, 6, 2)) == (-1/3, 1)

# (3x+1)^2 + 1 = 9x^2 + 6x + 2
@test unsafe_minimum(QuadraticPolynomial(9.0, 6.0, 2.0)) == 1.0

@test zero(QuadraticPolynomial{Float64}) == QuadraticPolynomial(0.0, 0.0, 0.0)

@test zero(QuadraticPolynomial(4.0, -2.0, 3.0)) == QuadraticPolynomial(0.0, 0.0, 0.0)

@test zero(QuadraticPolynomial(4, -2, 3)) == QuadraticPolynomial(0, 0, 0)
