using EllZeroTrendFiltering
using StaticArrays
using DelimitedFiles
using Plots

dpath = joinpath(dirname(@__FILE__), "data", "testo_data")
LH_27yo = readdlm(joinpath(dpath, "LH_27yo.csv"))
GnRH_27yo = readdlm(joinpath(dpath, "GnRH_27yo.csv"))
# TODO: also consider 40 yo

#b1 = 0.69 #Try different values for b1 and b2
b1 = 0.69
b2 = 0.014
Ac = @SMatrix[-b2 1; 0 -b1]
Bc = @SMatrix[0; 1/b1]
Cc = @SMatrix[1 0]

Ts = 10
A = exp(Ac*Ts)

# Obtaining discrete system, two approaches
#sysd, _ = c2d(ss(Matrix(Ac), Matrix(Bc), Matrix(Cc), 0), 10)
#sysd = ss(sysd.A, sysd.B, sysd.C, sysd.D, 1)
#sysd2 = ss(exp(Matrix(Ac)*10), Matrix(Bc), Matrix(Cc), 0)
#y, t, x = impulse(sysd, 20)
#y2, t2, x2 = impulse(sysd2, 20)

plot()

β2 = 2.0 # Basal secretion
g = LH_27yo[:, 2] .- β2
J_vec, U_vec, f_vec, y_vec = EllZeroTrendFiltering.fit_lti_constrained(exp(Ac*10), Cc, g, 20)
plot!(f_vec)



plot(0:10:1000, LH_27yo[:, 2], lw=2)
plot!(0:10:1000, y_vec[7] .+ β2, lw=2)


scatter(GnRH_27yo[:,1], GnRH_27yo[:,2], m=:o)
m = 8
scatter!(J_vec[m]*10, 0.35*U_vec[m])



# xo = [0, 0]
# t=0:10:1440

#t_true = [ 53 325 407 677 736 1007 1083 1352 1416 1686;  0.0015 0.0057 0.0063 0.0015 0.0059 0.0015 0.0075 0.0015 0.006 0.0065];
# u = zeros(1500)
# u[[53, 325, 407, 677, 736, 1007, 1083, 1352, 1416]] .= [0.0015, 0.0057, 0.0063, 0.0015, 0.0059, 0.0015, 0.0075, 0.0015, 0.006]
# y, t, x = lsim(sys, u, 1:length(u))
#
# dT = 5
# tred = 1:dT:length(y)
# g = y[tred]
# g = g .+ 0.1randn(N)
# N = length(g)
# sys_red = ss(exp(A*dT), B, C, 0, 1)

# plot(t, y)
# plot!(t2, y2)
