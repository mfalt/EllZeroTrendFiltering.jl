# PiecewiseLinearContinuousFitting
# Piecewise Linear Continouous Interpolation Constraints
# Sparse L2 Optimal Fitting subject Continuity Constraints
# Integrated Square

using Base.Test

include("../src/dev.jl")



using Polynomials
using IterTools

N = 120

t = linspace(0,2π,N)




#g = Poly([0,0,1,-1])
g = sin

@time ℓ = dev.compute_transition_costs(g, t)


K = 4
@time I, Y, f = dev.brute_force_optimization(ℓ, K-1)


Λ_0 = [dev.create_new_pwq(dev.minimize_wrt_x2(ℓ[i, N])) for i in 1:N-1]

@time Λ = dev.find_optimal_fit(Λ_0, ℓ, 7)

@time I2, y2, f2 = dev.recover_solution(Λ[3, 1], 1, N)
Y2, _ = dev.find_optimal_y_values(ℓ, I2)

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


#x = linspace(0,1,100)

#closefig()

using PyPlot
#close()
figure(1)
plot(t, g.(t), "g", t[I2], Y2, "bo-")

figure(2)
plot(1:length(t)-1, l')
