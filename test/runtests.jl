using EllZeroTrendFiltering
using Test, Random, LinearAlgebra
import Printf: @printf

tests = [
    "quadratic_polynomials.jl",
    "quadratic_forms.jl",
    "piecewise_quadratics.jl",
    "continuous_transition_costs.jl",
    "discrete_transition_costs.jl",
    "brute_force_search.jl",
    "pwq_dp_constrained.jl",
    "exact_tests.jl",
    "fit_pwl_constrained.jl",
    "fit_pwl_regularized.jl",
    "subsample_test.jl",
    "examples.jl"]

@testset "All Tests" begin
    @testset "Testfile: $test" for test in tests
        include(test)
    end
end
