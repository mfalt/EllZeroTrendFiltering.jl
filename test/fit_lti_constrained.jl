using EllZeroTrendFiltering
using StaticArrays
using ControlSystems


A = SMatrix{2,2,Float64,4}([0.9 0.1; 0 0.8])
B = SMatrix{2,1,Float64,2}([0; 0.2][:,:])
C = SMatrix{1,2,Float64,2}([1.0 0][:,:])

J_true = [4, 8, 10, 15]
U_true = [0.5, -0.5, 1, 2]

u_true = zeros(20)
u_true[J_true] .= U_true

sys = ss(Matrix{Float64}(A), [0; 1], Matrix{Float64}(C), 0, 1)
g0, _, _ = lsim(sys, u_true, 0:length(u_true)-1)

g = g0[:] + 0.2*randn(length(g0))

J_vec, U_vec, f_vec, y_vec = EllZeroTrendFiltering.fit_lti_constrained(A, C, g, 8) # FIXME: Doesn't work with 10

for m=1:4
    J_bf, U_bf, f_bf = brute_force_search(sys, g, m, initial_conditions=:zero)
    @test J_bf ≈ J_vec[m] atol=1e-10
    @test U_bf ≈ U_vec[m] atol=1e-10
    @test f_bf ≈ f_vec[m] atol=1e-10
end
