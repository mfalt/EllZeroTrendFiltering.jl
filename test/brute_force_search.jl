using Test
using EllZeroTrendFiltering

include("auxilliary_test_fcns.jl")

# Consider a subset of tests that are relatively small
for problem_fcn in ["discontinuous1",
                    "discontinuous2",
                    "super_exponential1",
                    "super_exponential3",
                    "white_noise",
                    "exponential"]

    include(joinpath(dirname(@__FILE__),"problems", problem_fcn * ".jl"))
    local g, I_sols  # TODO, can remove in julia 1.0?
    g, ζ_vec, I_sols, f_sols = @eval $(Symbol(problem_fcn))()

    M = length(I_sols)
    I_vec, _, f_vec = brute_force_multi(g, M)

    @testset "Data set: $problem_fcn, brute_force_search m=$m" for m in 1:length(f_sols)
        @test f_vec[m] ≈ f_sols[m]  atol=1e-8 # FIXME: The old value 1e-10 didn't work anymore
        if !isempty(I_sols[m])
            @test I_vec[m] == I_sols[m]
        end
    end
end
