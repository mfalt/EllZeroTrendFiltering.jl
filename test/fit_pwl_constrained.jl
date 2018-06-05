using Test
using EllZeroTrendFiltering


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

    include(joinpath(dirname(@__FILE__),"problems", problem_fcn * ".jl"))
    local g, I_sols  # TODO, can remove in julia 1.0?
    g, ζ_vec, I_sols, f_sols = @eval $(Symbol(problem_fcn))()

    M = length(I_sols)
    I_vec, _, f_vec = fit_pwl_constrained(g, M)

    @testset "Data set: $problem_fcn, constrained fit with m=$m" for m in 1:length(f_sols)
        @test f_vec[m] ≈ f_sols[m]  atol=1e-10
        if !isempty(I_sols[m])
            @test I_vec[m] == I_sols[m]
        end
    end

    # Iterate only over problems where solution exists
    @testset "Data set: $problem_fcn, subsampled, constrained fit with m=$m" for m in 1:length(f_sols)
        # Generate max(length(g)/5, m) random gridpoints in range 1:length(g)
        randI = shuffle(1:length(g))[1:max(floor(Int,length(g)/5), m+1)]
        # Merge with correct solution, and [1,length(g)] if solution doesn't exist
        t_grid = sort(union(randI, I_sols[m], [1, length(g)]))
        # Solve on subsampled grid,
        I_vec, _, f_vec = fit_pwl_constrained(g, m, t = t_grid)
        # Should be same if we had a solution
        if !isempty(I_sols[m])
            @test f_vec[m] ≈ f_sols[m]  atol=1e-10
            @test t_grid[I_vec[m]] == I_sols[m]
        else # Or higher cost if not
            @test f_vec[m] >= f_sols[m] || isapprox(f_vec[m], f_sols[m], atol=1e-10)
        end
    end
end
