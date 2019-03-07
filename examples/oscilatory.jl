using ControlSystems, EllZeroTrendFiltering
using Plots, Random, StaticArrays, DelimitedFiles
gr()

A = [0.993 0.05; -0.05 0.993]
B = [0;1]
C = [1 0]
sys = ss(A, B, C, 0, 1.0)

N = 1500
t = 0:1.0:N
u = zeros(N+1);

intime = [100, 250, 350, 450, 550, 650, 800, 900, 1100, 1300]
uvals = [3, -2, -1, 3, 1, 2, 3, -2, -3, 2]
u[intime] = uvals

y, t_out, x_out = lsim(sys, u, t)
y = y[:]

Random.seed!(2)
ynoise = y .+ 0.5.*randn(length(y))

@time J_vec, U_vec, f_vec, y_vec = EllZeroTrendFiltering.fit_lti_constrained(SMatrix{2,2}(A), SMatrix{1,2}(C), ynoise, 15)

lsimplot(sys, u, t, l=1, lab="true")
plot!(t, ynoise, c=:grey, opacity=0.5, lab="noisy")
plot!(t, y_vec[10], lab="10", l=(:dash,2))

writedlm("data_out_oscilatory/data.csv", [["gnoise"; ynoise] ["g"; y] ["y_l0"; y_vec[10][:]] ], ";")
#writedlm("data_out_oscilatory/impulse_data_l1.csv", [["T"; dSeq[1,:][:]] ["U"; dSeq[2,:][:]]], ";")
writedlm("data_out_oscilatory/impulse_data_l0.csv", [["T"; J_vec[10]]     ["U"; U_vec[10]]    ], ";")
writedlm("data_out_oscilatory/impulse_data_true.csv", [["T"; intime]     ["U"; uvals]    ], ";")


#plot!(t, y_vec[7], lab="7", l=(:dash,2))

plot(J_vec[10], U_vec[10], l=(:stem,:blue), m=(:circle, :blue));
plot!(intime, uvals, l=(:stem,:red), m=(:x, :red))

for i = 1:10
    uimp = zeros(length(t))
    uimp[intime[i]] = uvals[i]
    yimp = lsim(sys, uimp, t)[1]
    plot!(t, yimp[:], opacity=0.5, lab="");
end
plot!()
