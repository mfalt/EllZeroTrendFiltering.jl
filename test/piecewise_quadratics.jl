using Base.Test

include(joinpath(Pkg.dir("DynamicApproximations"),"src","dev.jl"))

pwq = dev.generate_PiecewiseQuadratic([([2.0, 2, 1], -Inf), ([2.0, -2, 1], 0.0)])

@test length(pwq) == 2
dev.add_quadratic2(pwq, dev.QuadraticPolynomial([1, 0, 1.0]))

@test length(pwq) == 3

#---
# Tests with (2x^2 + 1) and (x^2 + 2) which have intersections in -1, +1
# and 5/2*x(x-1) + 1 = 2.5x^2 - 2.5x + 1

p1 = dev.QuadraticPolynomial([2.0, 0, 1])
p2 = dev.QuadraticPolynomial([1.0, 0, 2])
p3 = dev.QuadraticPolynomial(2.5, -2.5, 1.0)

pwq1 = dev.create_new_pwq(p1)
dev.add_quadratic(pwq1, p2)
@test length(pwq1) == 3
dev.add_quadratic(pwq1, p3)
@test length(pwq1) == 4

pwq2 = dev.create_new_pwq(p2)
dev.add_quadratic(pwq2, p3)
@test length(pwq2) == 3
dev.add_quadratic(pwq2, p1)
@test length(pwq2) == 4

@test pwq1[1].p === pwq2[1].p == p1
@test pwq1[2].left_endpoint == pwq2[2].left_endpoint == -1
@test pwq1[2].p === pwq2[2].p == p2
@test pwq1[3].left_endpoint == pwq2[3].left_endpoint == 0
@test pwq1[3].p === pwq2[3].p

#---

# Given a matrix with coefficients for quadratic polynomials this function
# constructs the quadratic polynomials and inserts them into a piecewise
# quadratic object
function verify_piecewise_quadratic(c_mat)

    poly_list = [dev.QuadraticPolynomial(c_mat[k, :]) for k in 1:size(c_mat,1)]

    xgrid = linspace(-2,2)
    plot(show=true)
    for p in poly_list
        plot!(xgrid, p.(xgrid))
    end

    # Insert quadratics into piecewise quadfatic
    pwq = dev.create_new_pwq()
    for p in poly_list
        println("Hej")
        dev.add_quadratic2(pwq, p)
        println(pwq)
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

#---
c_mat1 = [
2.53106  1.59097   1.60448
2.51349  0.681205  1.60553
2.09364  0.714858  1.08607
]
verify_piecewise_quadratic(c_mat1)

#---
c_mat2 = [
2.01532  1.59269   1.65765
2.50071  1.56421   1.53899
2.02767  0.706143  1.129
]
verify_piecewise_quadratic(c_mat2)

#---

c_mat3 = [
#2.69548  1.50706   1.95186
#2.86858  0.990155  1.86534
2.42298  0.796953  1.02011
2.45055  1.32235   1.4894
2.1762   0.855014  1.0647
]
verify_piecewise_quadratic(c_mat3)
#---
# Random test case
N = 10
c_mat_rand = [2+rand(N) 2*rand(N) 1+rand(N)]

verify_piecewise_quadratic(c_mat_rand)
