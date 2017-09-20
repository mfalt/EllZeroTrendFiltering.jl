using DynamicApproximations
using Base.Test

tests = [
    "quadratic_polynomials.jl"
    "quadratic_forms.jl"
    "piecewise_quadratics.jl"
    "transition_costs.jl"
    "discrete_transition_costs.jl"
    "find_optimal_fit.jl"
    "test_discrete_fit.jl"]

@testset "All Tests" begin
    @testset "Testfile: $test" for test in tests
        include(test)
    end
end
