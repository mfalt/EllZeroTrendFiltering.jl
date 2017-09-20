using DynamicApproximations
using Base.Test

@testset "All Tests" begin
    include("quadratic_polynomials.jl")
    include("quadratic_forms.jl")
    include("piecewise_quadratics.jl")
    include("transition_costs.jl")
    include("discrete_transition_costs.jl")
    include("find_optimal_fit.jl")
    include("test_discrete_fit.jl")
end
