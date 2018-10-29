using Test
using EllZeroTrendFiltering: add_quadratic!

global const TEST_PRINT_LEVEL = 0

Random.seed!(27182)
#---

# Given a matrix with coefficients for quadratic polynomials this function
# constructs the quadratic polynomials and inserts them into a piecewise
# quadratic object
function verify_piecewise_quadratic(c_mat, show_plot=false)

    poly_list = [QuadraticPolynomial(c_mat[k, :]) for k in 1:size(c_mat,1)]

    if show_plot
        xgrid = linspace(-2,2)
        plot(show=true)
        for p in poly_list
            plot!(xgrid, p.(xgrid))
        end
    end

    # Insert quadratics into piecewise quadfatic
    pwq = create_new_pwq(Float64)
    for p in poly_list
        TEST_PRINT_LEVEL > 0 && println("Inserting: ", p)
        add_quadratic!(pwq, p)

        TEST_PRINT_LEVEL > 0 && println("After Insertion: ")
        TEST_PRINT_LEVEL > 0 && println(pwq)
    end

    # Check that the intersections between the quadratics are in increasing order
    TEST_PRINT_LEVEL > 0 && println("Checking that intersections points are increasing: ")
    for λ in pwq
        @test λ.left_endpoint < get_right_endpoint(λ)
    end

    # Check for continuity of the piecewise quadratic
    TEST_PRINT_LEVEL > 0 && println("Checking continuity: ")
    for λ in pwq
        x = get_right_endpoint(λ)
        if x == Inf
            break
        end

        TEST_PRINT_LEVEL > 0 && println("Continutiy diff: ", λ.p(x) - λ.next.p(x))
        @test λ.p(x) ≈ λ.next.p(x)
    end

    # Check that the piecewise quadratic is smaller than all the quadratics
    # that it was formed by
    x_grid = -5:0.1:5
    y_pwq = pwq(x_grid)
    # Test both versions of evaluating piecewise-quadratic
    @test y_pwq ≈ pwq.(x_grid)  rtol=sqrt(eps())

    # Test get_vals
    x, y, x_all, y_all = get_vals(pwq)
    @test pwq(x) ≈ y            rtol = sqrt(eps())

    for p in poly_list
        TEST_PRINT_LEVEL > 0 && @printf("suboptimality diff: %2.3f\n", maximum(y_pwq - p.(x_grid)))
        #println("   ",  x_grid[ (y_pwq - p.(x_grid)) .> 0])
        @test (y_pwq .<= p.(x_grid)) == trues(length(x_grid))
    end

end

pwq = generate_PiecewiseQuadratic(([2.0, 2, 1], -Inf), ([2.0, -2, 1], 0.0))

@test length(pwq) == 2
add_quadratic!(pwq, QuadraticPolynomial([1, 0, 1.0]))

#Test print string
str = "PiecewiseQuadratic{Float64} with 3 elements:\n[  -∞ , -0.50]\t  :   1.00*x^2 + 2.00*x + 2.00   \n[-0.50, 0.50]\t  :   1.00*x^2 + 0.00*x + 1.00   \n[0.50, ∞  ]\t  :   1.00*x^2 + -2.00*x + 2.00   \n"
str2 = sprint(print, pwq)
@test str == str2

@test length(pwq) == 3

#---
# Tests with (2x^2 + 1) and (x^2 + 2) which have intersections in -1, +1
# and 5/2*x(x-1) + 1 = 2.5x^2 - 2.5x + 1

p1 = QuadraticPolynomial([2.0, 0, 1])
p2 = QuadraticPolynomial([1.0, 0, 2])
p3 = QuadraticPolynomial(2.5, -2.5, 1.0)

pwq1 = create_new_pwq(p1)
add_quadratic!(pwq1, p2)
@test length(pwq1) == 3
add_quadratic!(pwq1, p3)
@test length(pwq1) == 4

pwq2 = create_new_pwq(p2)
add_quadratic!(pwq2, p3)
@test length(pwq2) == 3
add_quadratic!(pwq2, p1)
@test length(pwq2) == 4

@test pwq1[1].p === pwq2[1].p == p1
@test pwq1[2].left_endpoint == pwq2[2].left_endpoint == -1
@test pwq1[2].p === pwq2[2].p == p2
@test pwq1[3].left_endpoint == pwq2[3].left_endpoint == 0
@test pwq1[3].p === pwq2[3].p


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
    2.42298  0.796953  1.02011
    2.45055  1.32235   1.4894
    2.1762   0.855014  1.0647
    ]
verify_piecewise_quadratic(c_mat3)
#---
# Random test case
N = 10
c_mat_rand = [2 .+ rand(N) 2*rand(N) 1 .+ rand(N)]

verify_piecewise_quadratic(c_mat_rand)
