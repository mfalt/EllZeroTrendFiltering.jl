using DynamicApproximations
using Base.Test

tests = [
    "quadratic_polynomials.jl",
    "quadratic_forms.jl",
    "piecewise_quadratics.jl",
    "transition_costs.jl",
    "discrete_transition_costs.jl",
    "brute_force_search.jl",
    "pwq_dp_constrained.jl",
    "fit_pwl_constrained.jl",
    "fit_pwl_regularized.jl",
    "examples.jl"]

@testset "All Tests" begin
    @testset "Testfile: $test" for test in tests
        include(test)
    end
end
