using EllZeroTrendFiltering, Test

include(joinpath(dirname(@__FILE__),"auxilliary_test_fcns.jl"))

srand(1234)

# Brute force up to nbr segments:
M = 7
# Total number of points:
N = 40
# Number of internal gridpoints:
n = 12

#Function
g = cumsum(randn(N).^2)

t_grids = Vector{Int}[1:N, [[1;sort!(shuffle(2:N-1)[1:n]);N] for i = 1:10]...]
@testset "Testing subgrid, regularized, with endpoints, i=$i" for (i,t_grid) in enumerate(t_grids)
    local I_sols  # TODO, can remove in julia 1.0?
    I_sol, Y_sol, f_sol = brute_force_multi(g, M, t_grid, tol=1e-3)

    I, Y, f = fit_pwl_constrained(g, M, t=t_grid)
    I_t = getindex.([t_grid], I)
    @test I_t == I_sol
    @test Y ≈ Y_sol atol = 1e-10 # Crazy that this is possible
    @test f ≈ f_sol atol = 1e-10

    ζ_vec = 10 .^ range(log10(f_sol[2]), stop=log10(0.8*(max(f_sol[end-1] - f_sol[end], f_sol[end-1]))), length=3)

    for ζ in ζ_vec
        local I2, Y2  # TODO, can remove in julia 1.0?
        I2, Y2, f2 = fit_pwl_regularized(g, ζ, t=t_grid)
        I2_t = t_grid[I2]
        #Make sure regularized solution is in constrained solution
        idx = findall([I2_t] .== I_t)
        @test length(idx) == 1
        if length(idx) == 1
            @test Y2 ≈ Y_sol[idx[1]] atol = 1e-10
            @test f2 ≈ f_sol[idx[1]] atol = 1e-10
        end
    end
end

println("Testing sub-grid that is not covering the full range of input, expect warnings.")
# Test when endpoints are not included
t_grids = Vector{Int}[1:N, [sort!(shuffle(2:N-1)[1:n]) for i = 1:10]...]
@testset "Testing subgrid, constrained, without endpoints, i=$i" for (i, t_grid) in enumerate(t_grids)
    I_sol, Y_sol, f_sol = brute_force_multi(g, M-2, t_grid, tol=1e-3)

    I, Y, f = fit_pwl_constrained(g, M-2, t=t_grid)
    I_t = getindex.([t_grid], I)
    @test I_t == I_sol
    @test Y ≈ Y_sol atol = 1e-10 # Crazy that this is possible
    @test f ≈ f_sol atol = 1e-10

    ζ_vec = 10 .^ range(log10(f_sol[2]), stop=log10(0.8*(max(f_sol[end-1] - f_sol[end], f_sol[end-1]))), length=3)

    for ζ in ζ_vec
        local I2, Y2  # TODO, can remove in julia 1.0?
        I2, Y2, f2 = fit_pwl_regularized(g, ζ, t=t_grid)
        I2_t = t_grid[I2]
        #Make sure regularized solution is in constrained solution
        idx = findall([I2_t] .== I_t)
        @test length(idx) == 1
        if length(idx) == 1
            @test Y2 ≈ Y_sol[idx[1]] atol = 1e-10
            @test f2 ≈ f_sol[idx[1]] atol = 1e-10
        end
    end
end
