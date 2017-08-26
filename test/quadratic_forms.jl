using Base.Test

using StaticArrays

include("../src/dev.jl")


# Addition of quaratic forms
Q = dev.QuadraticForm2(@SMatrix([2.0 0; 0 1]), @SVector([0.0,0]), 1.0) +
dev.QuadraticForm2(@SMatrix([2.0 1; 1 2]), @SVector([1.0,0]), 0.0)

@test Q == dev.QuadraticForm2(@SMatrix([4.0 1; 1 3]), @SVector([1.0,0]), 1.0)



# Maybe somewhat surprising that the tests pass considering equalities between
# Floats, or is the aritmhetic exact

# x^2 + y^2 + 1
Q1 = dev.QuadraticForm2(@SMatrix([1 0; 0 1]), @SVector([0, 0]), 1)
@test dev.minimize_wrt_x2(Q1) == dev.QuadraticPolynomial(1,0,1)
@test Q1(1, 2) == 6.0


# (2x - y + 1)^2 + x^2 - 2x + 1
# = 5x^2 + y^2 + 2 + 2x - 2y - 4xy
Q2 = dev.QuadraticForm2(@SMatrix([5 -2; -2 1]) ,
@SVector([2, -2]),
2)
@test dev.minimize_wrt_x2(Q2) == dev.QuadraticPolynomial(1,-2,1)
@test Q2(1, 2) - 1.0 < 1e-14


# x^2 + 2x + 2
Q3 = dev.QuadraticForm2(@SMatrix([1 0; 0 0]),
@SVector([2, 0]),
2)
@test dev.minimize_wrt_x2(Q3) == dev.QuadraticPolynomial(1,2,2)


# y^2 + 2y + 2
Q4 = dev.QuadraticForm2(@SMatrix([0 0; 0 1]),
@SVector([0, 2]),
2)

@test dev.minimize_wrt_x2(Q4) == dev.QuadraticPolynomial(0,0,1)
