using Test
using EllZeroTrendFiltering
using Plots
plotly()

include("write_problem_to_file.jl")
include("auxilliary_test_fcns.jl")

problem_data = Vector{Tuple{String,Any,Integer}}(undef, 10)

problem_data[1] = ("discontinuous1",
    [0.2*linear_trend(10) + circle_segment(10);
     -6 + 8*linear_trend(5);
     -6 + 8*linear_trend(5)], 9)
# For m=8, there are multiple solutions that are within numerical precision

problem_data[2] = ("discontinuous2",
    [5*circle_segment(7);
    -5 + 4*linear_trend(5);
    -5 + 4*linear_trend(5);
     5*circle_segment(7)], 8)


# A234588
problem_data[3] = ("OEIS1",
     [1.0, 2, 4, 5, 9, 14, 8, 10, 18, 27, 17, 6, 11, 15, 29, 44, 19, 36, 54, 35,
     21, 42, 23, 46, 25, 50, 28, 55, 83, 112, 31, 62, 33, 66, 37, 72, 108, 71,
     39, 78, 41, 82, 124, 45, 89, 134, 88, 48, 96, 51, 101, 152, 53, 106, 160],
     5)

problem_data[4] = ("white_noise",
    randn(srand(1), 15),
    14)

problem_data[5] = ("exponential",
    exp.(0:0.25:4),
    16)

problem_data[6] = ("super_exponential1",
    exp.((0:0.25:3).^2),
    12)

problem_data[7] = ("super_exponential2",
    [exp.((0:0.25:2).^2) ; exp((2:-0.25:0.5).^2)],
    15)

problem_data[8] = ("super_exponential3",
    [exp.((0:0.25:2).^2) ; -exp((2:-0.25:0.5).^2)],
    15)

problem_data[9] = ("super_exponential4",
    [exp((2:-0.25:0.5).^2) ; exp.((0:0.25:2).^2)],
    15)

problem_data[10] = ("super_exponential5",
    [exp((2:-0.25:0).^2) ; 2-exp.((0:0.25:2).^2)],
    17)



for (problem_name, g, M) in problem_data[10:end]


    @time I_vec, Y_vec, f_vec = brute_force_multi(g, M)

    # It should not be beneficial to use more segments than has been computed ...
    ζ_vec = 10 .^ range(log10(f_vec[2]), stop=log10(0.8*(max(f_vec[end-1] - f_vec[end], f_vec[end-1]))), length=10)

    write_problem_to_file(problem_name, g, ζ_vec, I_vec, f_vec)

    plot(g, m="o")
    for k=1:length(I_vec)
        plot(I_vec[k]-1, Y_vec[k])
    end

end

## provide some cheap options to compute costs
# using Interpolations
#
# ℓ = compute_discrete_transition_costs(g)
# V_N = QuadraticPolynomial(1.0, -2*g[end], g[end]^2)
#
# Y, f = EllZeroTrendFiltering.find_optimal_y_values(ℓ, V_N, I)
#
# y = interpolate((I,), Y, Gridded(Linear()))(1:length(g))
# cost2 = sum((y[1:end]-g[1:end]).^2) # Note: cost at i=N should not be included
