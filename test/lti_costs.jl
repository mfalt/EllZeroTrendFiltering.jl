### Double integrator system..
A = SMatrix{2,2,Float64,4}([1 1; 0 1])
B = SMatrix{2,1,Float64,2}([0; 1])
C = SMatrix{1,2,Float64,2}([1 0])

# Constant g with all ones
g = ones(5)
l, χ, V_N = compute_problem_data_lti(g, A, C)

@test χ[1,2] == [1 1]
@test χ[1,3] == [1 2]
@test χ[1,4] == [1 3]


@test l[1,5]([0,0]) == 4
@test l[1,5]([0,1]) == 1^2 + 0^2 + 1^2 + 2^2
@test l[1,5]([1,0]) == 0
@test l[1,5]([0.5,0]) == 1

@test l[1,5]([1,1]) == 0^2 + 1^2 + 2^2 + 3^2
@test l[3,5]([1,1]) == 0^2 + 1^2


# Straight line
g = Vector{Float64}(0:4)
l, _ = compute_problem_data_lti(g, A, C)

# For double integrator system..
@test l[1,5]([0,0]) == 0^2 + 1^2 + 2^2 + 3^2
@test l[1,5]([1,0]) == 1^2 + 0^2 + 1^2 + 2^2
@test l[3,5]([2,0]) == 0^2 + 1^2
@test l[1,5]([0,1]) == 0
@test l[1,5]([1,1]) == 4



### General serond order system
A = SMatrix{2,2,Float64,4}([0.9 0.1; 0 0.8])
B = SMatrix{2,1,Float64,2}([0; 0.2][:,:])

# g is all zeros
g = zeros(10)
l, χ_vec = compute_problem_data_lti(g, A, C)
V_N = QuadraticPolynomial(1.0, -2*g[end], g[end]^2)

@test χ_vec[1,2] == [0.9 0.1]
@test χ_vec[1,3] == [0.9^2 0.1*(0.9+0.8)]
@test χ_vec[1,4] ≈ [0.9^3 0.1*(0.9^2+0.9*0.8+0.8^2)]
@test χ_vec[1,5] ≈ [0.9^4 0.1*(0.9^3+0.9^2*0.8+0.9*0.8^2+0.8^3)]

@test l[1,3]([0.0,1]) == 0.1^2

sys = ss(Matrix{Float64}(A), Matrix{Float64}(B), Matrix{Float64}(I, 2, 2), 0.0, 1)
y, t, x = lsim(sys, zeros(10), 0:9, [0.0, 1])
x = x'

j = 5
@test sum(abs.(g[1:j-1] - x[1,1:j-1]).^2) == l[1,j]([0.0,1])


# g is all zeros but multiple impulses
x0 = [0.0, 1]
u = zeros(10)
u[1] = 1.0
u[5] = -1.0
y, t, x = lsim(sys, zeros(10), 0:9, [0.0, 1])
x = x'

@test l[1,2](x[:,1]) + l[2,6](x[:,2]) + l[6,10](x[:,6]) ≈
    sum(abs.(g[1:9] - x[1,1:9]).^2)

@test l[1,2](x[:,1]) + l[2,6](x[:,2]) + l[6,10](x[:,6]) + V_N(x[1,10]) ≈
    sum(abs.(g[1:10] - x[1,1:10]).^2)



# g is more general
g = [1.0, 0.0, 1.0, 2.0, 1.5, 1.5, 1.0, 1.0, 0.5, 0.0]
l, _ = compute_problem_data_lti(g, A, C)
V_N = QuadraticPolynomial(1.0, -2*g[end], g[end]^2)

x0 = [0.0, 1]
u = zeros(10)
u[1] = 1.0
u[5] = -1.0

sys = ss(Matrix{Float64}(A), Matrix{Float64}(B), Matrix{Float64}(I, 2, 2), 0.0, 1)
y, t, x = lsim(sys, zeros(10), 0:9, [0.0, 1])
x = x'
#x = lsim(A, B, u, x0)

@test l[1,2](x[:,1]) + l[2,6](x[:,2]) + l[6,10](x[:,6]) ≈
    sum(abs.(g[1:9] - x[1,1:9]).^2)

@test l[1,2](x[:,1]) + l[2,6](x[:,2]) + l[6,10](x[:,6]) + V_N(x[1,10]) ≈
    sum(abs.(g[1:10] - x[1,1:10]).^2)
