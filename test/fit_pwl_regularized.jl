using Base.Test
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

    include(joinpath(Pkg.dir("EllZeroTrendFiltering"),"test","problems", problem_fcn * ".jl"))

    g, ζ_vec, I_sols, f_sols = @eval $(Symbol(problem_fcn))()

    ζ = 0.01

    @testset "Data set: $problem_fcn, regularization with ζ=$ζ" for ζ in ζ_vec

        I_reg, _, f_reg = fit_pwl_regularized(g, ζ)
        # (the cost f_reg includes penality mζ)

        # Use costs in the solution file to find out how many segments the solution should contain
        m_expected = indmin([f_sols[m] + m*ζ for m=1:length(f_sols)])

        @test m_expected == length(I_reg) - 1
        if !isempty(I_sols[m_expected]) # Only if there is a unique solution
            @test I_reg == I_sols[m_expected]
        end
        @test f_reg ≈ f_sols[m_expected]   atol=1e-10  # + ζ*m_expected


    end

    # Iterate only over problems where solution exists
    @testset "Data set: $problem_fcn, subsampled, regularization with ζ=$ζ" for ζ in ζ_vec
        # Use costs in the solution file to find out how many segments the solution should contain
        m_expected = indmin([f_sols[m] + m*ζ for m=1:length(f_sols)])

        # Generate max(length(g)/5, m) random gridpoints in range 1:length(g)
        randI = shuffle(1:length(g))[1:max(floor(Int,length(g)/5), m_expected+1)]
        # Merge with correct solution, and [1,length(g)] if solution doesn't exist
        t_grid = sort(union(randI, I_sols[m_expected], [1, length(g)]))
        # Solve on subsampled grid
        I_reg, _, f_reg = fit_pwl_regularized(g, ζ; t = t_grid)

        # Should be same if we had a solution
        if !isempty(I_sols[m_expected])
            @test m_expected == length(I_reg) - 1
            @test f_reg ≈ f_sols[m_expected]  atol=1e-10
            @test t_grid[I_reg] == I_sols[m_expected]
        else
            @test f_reg >= f_sols[m_expected] || isapprox(f_reg, f_sols[m_expected], atol=1e-10)
        end
    end
end
