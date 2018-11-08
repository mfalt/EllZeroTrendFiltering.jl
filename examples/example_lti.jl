##
using Revise
using EllZeroTrendFiltering
using ControlSystems
using StaticArrays
using LinearAlgebra



A = SMatrix{2,2,Float64,4}([0.9 0.1; 0 0.8])
B = SMatrix{2,1,Float64,2}([0; 0.2][:,:])
C = SMatrix{1,2,Float64,2}([1.0 0][:,:])


#u_true[3] = 2
#u_true[15] = -1
J_true = [4, 8, 22, 40]
U_true = [0.5, -0.5, 1, 2]

u_true = zeros(59)
u_true[J_true] .= U_true
u_true = u_true[1:30]

sys = ss(Matrix{Float64}(A), [0; 1], Matrix{Float64}(C), 0, 1)
g0, t, x = lsim(sys, u_true, 0:length(u_true)-1)

g0 = 2g0[:]
g = g0 + 0.05*randn(length(g0))

N = length(g)

l, χ, V_N = compute_problem_data_lti(g, A, C)

Λ = construct_value_fcn_constrained(l, χ, V_N, 12, 1000)

m = 10

# Zero initial conditions
#J_set0, f = recover_optimal_index_set_zero_ic(l, Λ, m)

J, f = EllZeroTrendFiltering.recover_optimal_index_set(Λ[8, :], l, χ, :zero)

J
U_vec, y_vec = EllZeroTrendFiltering.post_process_lti(A, C, g, [J])
y = y_vec[1]; U = U_vec[1]

norm(y - g)^2


J_vec, U_vec, f_vec, y_vec = EllZeroTrendFiltering.fit_lti_constrained(A, C, g, 3)

for m=1:4
    J_opt, U_opt, f_opt = brute_force_search(sys, g, m, initial_conditions=:zero)
    println(f_opt)
end





# Positiva impulser


using Plots


#gui()

plot(g, m=:o,l=nothing)
plot!(y_vec[4], m=:o, color="blue", l=nothing)


plot(f_vec)

#plot!(Y_opt, m=:o, color="red", l=nothing)


# Lsim discrete systems, shouldn't have to specify time
# State space of static arrays


# Kan droppa en del typparametrar i SVector, i.e. only SMatrix{2,2,T}, SMatrix{2,2}?

##
# Testa med 0 intialvillkor,
#
# Syntetiskt genererad data, data + brus



f_opt^2
