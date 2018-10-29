# PiecewiseLinearContinuousFitting
# Piecewise Linear Continouous Interpolation Constraints
# Sparse L2 Optimal Fitting subject Continuity Constraints
# Integrated Square
using EllZeroTrendFiltering: minimize_wrt_x2

using Polynomials
using IterTools

N = 240

t = range(0, stop=4π, length=50)




#g = Poly([0,0,1,-1])
ran = 0.1*randn(N)
z(t) = ran[floor(Int64, t/4pi*(N-1))+1]
g = sin
g2(t) = g(t) + z(t)

@time ℓ = compute_transition_costs(g, t);


K = 7
#@time I, Y, f = brute_force_search(ℓ, K-1);


Λ_0 = [create_new_pwq(minimize_wrt_x2(ℓ[i, N])) for i in 1:N-1];

@time Λ = pwq_dp_constrained(Λ_0, ℓ, 7);
tot

using Plots
plot(t, g.(t), lab="g(t)")
for k = 7:7
    @time I2, y2, f2 = recover_solution(Λ[k, 1], 1, N)
    Y2, _ = find_optimal_y_values(ℓ, I2)

    println("Comparison: ", sqrt(f), " --- ", sqrt(f2))


    #find_optimal_y_values


    l = zeros(size(Λ))
    for i=1:size(l,1)
        for j=1:size(l,2)
            if isassigned(Λ,i,j)
                l[i,j] = length(Λ[i, j])
            end
        end
    end


    #x = range(0, stop=1, length=100)

    #closefig()
    plot!(t[I2],Y2, m=:circle, lab="k=$k")
end
plot!()
#using PyPlot
# #close()
# figure(1)
# plot(t, g.(t), "g", t[I2], Y2, "bo-")
#
# figure(2)
# plot(1:length(t)-1, l')
