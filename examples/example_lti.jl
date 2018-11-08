##
using Revise
using EllZeroTrendFiltering
using ControlSystems
using StaticArrays
using LinearAlgebra

function compute_response(A, C, N, J, U; x0=[0; 0])
    u = zeros(N)
    u[J] .= U
    sys = ss(Matrix{Float64}(A), [0; 1], Matrix{Float64}(C), 0, 1)
    return lsim(sys, u, 1:N)[1]
end

A = SMatrix{2,2,Float64,4}([1.0 1.0; 0 1.0])
B = SMatrix{2,1,Float64,2}([0; 1.0][:,:])
C = SMatrix{1,2,Float64,2}([1.0 0][:,:])

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


N = length(u0)
sys = ss(Matrix{Float64}(A), [0; 1], Matrix{Float64}(C), 0, 1)
g0, t, x = lsim(sys, u_true, 0:length(u_true)-1)

g0 = 2g0[:]
g = g0 + 0.05*randn(length(g0))

N = length(g)

l, χ, V_N = compute_problem_data_lti(g, A, C)
#l, χ, V_N = compute_problem_data_pwl(g, 1:N)
Λ = construct_value_fcn_constrained(l, χ, V_N, 8, 1000)




m = 4

find_minimum(Λ[4, 1])
J_set0, f = recover_optimal_index_set_zero_ic(l, Λ, m)

# Zero initial conditions
c = 0; X = generate_markov_matrix(sys, N, free_intial_conditions=false)
Xs = X[:, [1:c; J_set0 .- 1]]
U_opt = Xs\g
norm(Xs*U_opt - g)^2

y = compute_response(A, C, length(g), J_set0 .- 1, U_opt)[:]
norm(y - g)^2


J_opt, U_opt, f_opt = brute_force_search(X, g, m, 0)


# Free initial conditions
J_set_free, y, f = recover_optimal_index_set(Λ[m,1])

c = 2; X = generate_markov_matrix(sys, N, free_intial_conditions=true)
Xs = X[:, [1:c; J_set_free .- 1 .+ c]]
U_opt = Xs\g
norm(Xs*U_opt - g)^2
y = compute_response(A, C, length(g), J_set0 .- 1, U_opt)
norm(Xs*U_opt - g)^2

J_opt, U_opt, f_opt = brute_force_search(X, g, m-1, 2)


# Positiva impulser


using Plots


#gui()

plot(g, m=:o,l=nothing)
#plot!(g0, m=:o,c="green",l=nothing)
#plot!(A_markov[:,J]*th, m=:o,l=nothing)
plot!(y, m=:o, color="blue", l=nothing)
plot!(Y_opt, m=:o, color="red", l=nothing)


# Lsim discrete systems, shouldn't have to specify time
# State space of static arrays


# Kan droppa en del typparametrar i SVector, i.e. only SMatrix{2,2,T}, SMatrix{2,2}?

##
# Testa med 0 intialvillkor,
#
# Syntetiskt genererad data, data + brus



f_opt^2
