p0 = QuadraticPolynomial(0.0, 0.0, 0.0)
@test (p0 ∘ @SMatrix[1 0]) == QuadraticForm{Float64}([0.0 0.0; 0.0 0.0], [0.0, 0.0], 0.0)
@test (p0 ∘ @SMatrix[0 1]) == QuadraticForm{Float64}([0.0 0.0; 0.0 0.0], [0.0, 0.0], 0.0)

p1 = QuadraticPolynomial(0.0, 0.0, 1.0)
@test (p1 ∘ @SMatrix[1 0]) == QuadraticForm{Float64}([0.0 0.0; 0.0 0.0], [0.0, 0.0], 1.0)
@test (p1 ∘ @SMatrix[0 1]) == QuadraticForm{Float64}([0.0 0.0; 0.0 0.0], [0.0, 0.0], 1.0)

p2 = QuadraticPolynomial(0.0, -1.0, 1.0)
@test (p2 ∘ @SMatrix[1 0]) == QuadraticForm{Float64}([0.0 0.0; 0.0 0.0], [-1.0, 0.0], 1.0)
@test (p2 ∘ @SMatrix[0 1]) == QuadraticForm{Float64}([0.0 0.0; 0.0 0.0], [0.0, -1.0], 1.0)

p3 = QuadraticPolynomial(2.0, 0.0, 1.0)
@test (p3 ∘ @SMatrix[1 0]) == QuadraticForm{Float64}([2.0 0.0; 0.0 0.0], [0.0, 0.0], 1.0)
@test (p3 ∘ @SMatrix[0 1]) == QuadraticForm{Float64}([0.0 0.0; 0.0 2.0], [0.0, 0.0], 1.0)
@test (p3 ∘ @SMatrix[2 1]) == QuadraticForm{Float64}([8.0 4.0; 4.0 2.0], [0.0, 0.0], 1.0)
