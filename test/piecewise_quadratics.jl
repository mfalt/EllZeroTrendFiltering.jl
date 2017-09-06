using Base.Test

include(joinpath(Pkg.dir("DynamicApproximations"),"src","dev.jl"))

pwq = dev.generate_PiecewiseQuadratic([([2.0, 2, 1], -Inf), ([2.0, -2, 1], 0.0)])

@test length(pwq) == 2
dev.add_quadratic2(pwq, dev.QuadraticPolynomial([1, 0, 1.0]))

@test length(pwq) == 3

#---
# Tests with (2x^2 + 1) and (x^2 + 2) which have intersections in -1, +1

p1 = dev.QuadraticPolynomial([2.0, 0, 1])
p2 = dev.QuadraticPolynomial([1.0, 0, 2])

pwq1 = dev.create_new_pwq(p1)
dev.add_quadratic2(pwq1, p2)

println(pwq1)



pwq2 = dev.create_new_pwq(p2)
dev.add_quadratic(pwq2, p1)

println(pwq1)

#---

# Given a matrix with coefficients for quadratic polynomials this function
# constructs the quadratic polynomials and inserts them into a piecewise
# quadratic object
function verify_piecewise_quadratic(c_mat)

    poly_list = [dev.QuadraticPolynomial(c_mat[k, :]) for k in 1:size(c_mat,1)]

    # Insert quadratics into piecewise quadfatic
    pwq = dev.create_new_pwq()
    for p in poly_list
        dev.add_quadratic2(pwq, p)
    end

    # Check that the intersections between the quadratics are in increasing order
    for λ in pwq
        @test λ.left_endpoint < dev.get_right_endpoint(λ)
    end

    # Check for continuity of the piecewise quadratic
    for λ in pwq
        if isnull(λ.next)
            break
        end
        x = get(λ.next).left_endpoint
        println("Continutiy diff: ", λ.p(x) - get(λ.next).p(x))
        @test λ.p(x) ≈ get(λ.next).p(x)
    end

    # Check that the piecewise quadratic is smaller than all the quadratics
    # that it was formed by
    x_grid = -5:0.1:5
    y_pwq = dev.evalPwq(pwq, x_grid)

    for p in poly_list
        @printf("suboptimality diff: %2.3f\n", maximum(y_pwq - p.(x_grid)))
        #println("   ",  x_grid[ (y_pwq - p.(x_grid)) .> 0])
        @test (y_pwq .<= p.(x_grid)) == trues(length(x_grid))
    end

end


c_mat1 = [
2.53106  1.59097   1.60448
2.51349  0.681205  1.60553
2.09364  0.714858  1.08607
]
verify_piecewise_quadratic(c_mat1)

c_mat2 = [
2.01532  1.59269   1.65765
2.50071  1.56421   1.53899
2.02767  0.706143  1.129
]
verify_piecewise_quadratic(c_mat2)


# Random test case
N = 10
c_mat_rand = [2+rand(N) 2*rand(N) 1+rand(N)]

verify_piecewise_quadratic(c_mat_rand)
