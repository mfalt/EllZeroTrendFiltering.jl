using Base.Test, StaticArrays, DynamicApproximations
using DynamicApproximations.minimize_wrt_x2

# Addition of quaratic forms
Q = QuadraticForm(@SMatrix([2.0 0; 0 1]), @SVector([0.0,0]), 1.0) +
    QuadraticForm(@SMatrix([2.0 1; 1 2]), @SVector([1.0,0]), 0.0)
@testset "Quadratic Forms" begin
    @test Q == QuadraticForm(@SMatrix([4.0 1; 1 3]), @SVector([1.0,0]), 1.0)



    # Maybe somewhat surprising that the tests pass considering equalities between
    # Floats, or is the aritmhetic exact

    # x^2 + y^2 + 1
    Q1 = QuadraticForm(@SMatrix([1.0 0; 0 0]), @SVector([0.0, 0]), 0.5)
    p1 = QuadraticPolynomial(1.0, 0.0, 0.5)
    @test minimize_wrt_x2(Q1, p1) == QuadraticPolynomial(1,0,1)
    @test Q1(1, 2) == 1.5


    # (2x - y + 1)^2 + x^2 - 2x + 4
    # = 5x^2 + y^2 + 2 + 2x - 2y - 4xy
    Q2 = QuadraticForm(@SMatrix([5.0 -2; -2 0.5]), @SVector([2.0, -2]), 1.0)
    p2 = QuadraticPolynomial(0.5, 0.0, 5.0)
    @test minimize_wrt_x2(Q2, p2) == QuadraticPolynomial(1,-2,5)

    # x^2 + 2x + 2
    Q3 = QuadraticForm(@SMatrix([1.0 0; 0 0]), @SVector([2.0, 0]), 2.0)
    p3 = QuadraticPolynomial(0.0, 0.0, 0.0)
    @test minimize_wrt_x2(Q3, p3) == QuadraticPolynomial(1.0,2.0,2.0)


    # y^2 + 2y + 2
    Q4 = QuadraticForm(@SMatrix([0.0 0; 0 1]), @SVector([0.0, 2]),2.0)
    p4 = QuadraticPolynomial(0.0, 0.0, 0.0)
    @test minimize_wrt_x2(Q4, p4) == QuadraticPolynomial(0,0,1)
end
