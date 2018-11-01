using Plots
using Convex
using Mosek
using EllZeroTrendFiltering

using Interpolations


g = snp500_data()[1:2000]

N = length(g)

# Transition costs (for recvoering the y-values)
l = compute_discrete_transition_costs(g)
V_N = QuadraticPolynomial(1.0, -2*g[end], g[end]^2)

M = 50

#f = open("data.jld", "w+"); serialize(f, Λ); close(f)
#@time Λ =  construct_value_fcn_constrained(l, V_N, M)

sols_opt = []
for m=1:M
    I_opt, Y_opt, f_opt = recover_solution(Λ[m, 1], l, V_N)
    push!(sols_opt, (I_opt, Y_opt, f_opt))
end

# Iterative re-weighting according to Kim 2009, et. al
function iterative_l1_reweighting(g, max_cost, ϵ)
    N = length(g)
    H =  spdiagm((ones(N-2), -2*ones(N-2), ones(N-2)), (0,1,2))

    W = spdiagm(ones(N-2))

    x_sol = zeros(size(N))
    f_sol = NaN

    for iter=1:10
        x = Variable(N)
        problem = minimize(norm((W*H)*x, 1))
        problem.constraints += sumsquares(x - g) < max_cost

        solve!(problem, MosekSolver(), verbose=false)
        println(problem.status)

        x_sol = evaluate(x)[:]
        f_sol = sum((x_sol - g).^2)

        W = spdiagm(ones(N-2) ./ (ϵ .+ H*x_sol))
    end

    I_sol = [1; find( abs.(H*x_sol) .> 0.01ϵ); N]

    return (x_sol, I_sol, f_sol)
end


plot(g)
opt_f = []


ϵ_vec = [1e-5, 1e-6]
max_cost_vec = 10 .^ range(log10(0.1), stop=log10(8), length=2)


sols_l1_mat = []
for (i, max_cost)=enumerate(max_cost_vec)
    I = []
    Y = []
    f = []
    for (j, ϵ) = enumerate(ϵ_vec)
        _, I_new, f_new = iterative_l1_reweighting(g, max_cost, ϵ)

        if !isempty(I)
            if I != I_new
                println((i, j))
                println(I)
                println(I_new)
                println(f, ": ", f_new)
                _ = readline()
            end
        end

        I = I_new

        print(I)

        # Find optimal y values and optimal cost
        Y, f = EllZeroTrendFiltering.find_optimal_y_values(l, V_N, I)

    end
    push!(sols_l1_mat, (I, Y, f))
end


I_sols_opt_vec, Y_sols_opt_vec, f_sols_opt_vec = collect(zip(sols_opt...))
f_sols_opt_vec = [f_sols_opt_vec...]

I_sols_l1_vec, Y_sols_l1_vec, f_sols_l1_vec = collect(zip(sols_l1_mat...))
f_sols_l1_vec = [f_sols_l1_vec...]

m_sols_l1_vec = [length.(I_sols_l1_vec)...]
m_sols_opt_vec = [length.(I_sols_opt_vec)...]

plot(m_sols_l1_vec, f_sols_l1_vec);
plot!(m_sols_opt_vec, f_sols_opt_vec, xlims=(0,30), ylims=(0,10));
gui()



plot(g);
#plot!(x_sol, lw=2.5)
plot!(I_vec[m], Y_vec[m], lw=2.5)

plot(H*x_sol);
plot!(I_vec[m][2:end], 0.1*diff(Y_vec[m]) ./ diff(I_vec[m]), m=:o, lw=0)


writedlm("cost_vs_m_opt.csv", [["m"; m_sols_opt_vec] ["cost"; f_sols_opt_vec]], ";")
writedlm("cost_vs_m_l1.csv", [["m"; m_sols_l1_vec] ["cost"; f_sols_l1_vec]] , ";")



writedlm(".csv", [["I"; I_sols_opt_vec[10]] ["Y"; Y_sols_opt_vec[10]]] , ";")



# Computation times for ℓ_0 trend filtering
@time fit_pwl_regularized(g, 0.2)
@time fit_pwl_regularized(g, 0.1)
@time fit_pwl_regularized(g, 0.05)
@time fit_pwl_regularized(g, 0.01)
@time fit_pwl_regularized(g, 0.2)
