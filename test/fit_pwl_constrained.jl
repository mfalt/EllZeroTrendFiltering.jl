using Base.Test
using DynamicApproximations


for problem_fcn in ["discontinuous1",
                    "discontinuous2",
                    "OEIS1",
                    "super_exponential1",
                    "super_exponential2",
                    "super_exponential3",
                    "super_exponential4",
                    "super_exponential5",
                    "square_wave",
                    "white_noise",
                    "exponential",
                    "straight_line"]

    include(joinpath(Pkg.dir("DynamicApproximations"),"test","problems", problem_fcn * ".jl"))

    g, ζ_vec, I_sols, f_sols = @eval $(Symbol(problem_fcn))()

    M = length(I_sols)
    I_vec, _, f_vec = fit_pwl_constrained(g, M)

    @testset "Data set: $problem_fcn, constrained fit with m=$m" for m in 1:length(f_sols)
        @test f_vec[m] ≈ f_sols[m]  atol=1e-10
        if !isempty(I_sols[m])
            @test I_vec[m] == I_sols[m]
        end
    end
end
