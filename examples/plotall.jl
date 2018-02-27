
using Base.Test

@everywhere include(joinpath(Pkg.dir("EllZeroTrendFiltering"),"src","jl"))

@everywhere using Polynomials, IterTools, Plots


function plotall(Λ,ℓ,t,g, both=true)
    M = size(Λ,1)
    N = length(t)
    γ = both ? 2 : 1
    p = plot(layout=(M*γ,N-1), size=(1400,1000))
    k = 1
    for m = M:-1:1
        for n = 1:N-1
            idx = (n-1)*M+m #Linear index in Λ
            if isassigned(Λ,idx)
                x, y, x_all, y_all = get_vals(Λ[m,n])
                for i in eachindex(x_all)
                    plot!(p, y_all[i][[1,end]], x_all[i][[1,end]], l = 0, m=(1,:cross,:orange), subplot=k, lab="")
                    plot!(p, y_all[i], x_all[i], l=(1,:blue), lab="", subplot=k)
                end
                I, yI, fI = recover_solution(Λ[m,n], n, N)
                ystr = @sprintf("%3.2f", yI)
                fstr = @sprintf("%3.2f", fI)
                plot!(p, y, x, l=:red, subplot=k, lab="$I, $ystr, $fstr", ylims=(-4,4), xlims=(0,6))
                if both && length(I) > 1 && fI != NaN
                    Y, _ = find_optimal_y_values(ℓ, I)
                    tcont = linspace(t[1],t[end], 100*length(t))
                    plot!(p, tcont/t[end]*(N-1)+1, g.(tcont), l=:blue, lab="", subplot=k+N-1)
                    plot!(p, I, Y, m=:circle, l=:red, lab="", subplot=k+N-1)
                end
            end
            k += 1
        end
        k += both ? N-1 : 0
    end
    return p
end

##
srand(9)
N = 20
ran = cumsum(0.1*randn(10*N))

t = linspace(0,1,N)

g6(t) = ran[floor(Int64, t*10(N-1)+1)]
ℓ = compute_transition_costs(g6, t);

Λ_0 = [create_new_pwq(minimize_wrt_x2(ℓ[i, N], QuadraticPolynomial{Float64}(0.,0.,0.))) for i in 1:N-1];
Λ, t2, _, _, _ = @timed pwq_dp_constrained(Λ_0, ℓ, 19);
##
#[(isassigned(Λ,i) ? length(Λ[i]) : 0) for i in eachindex(Λ)]
lengths = reshape([(isassigned(Λ,i) ? length(Λ[i]) : 0) for i in eachindex(Λ)], size(Λ))
plotly()
heatmap(lengths, show=true)
# p = plotall(Λ,ℓ,t,g, true);
# @time I2, y2, f2 = recover_solution(Λ[19, 1], 1, N)
# Y2, _ = find_optimal_y_values(ℓ, I2)
# plot(linspace(1,N,20*N), g5.(linspace(1,N,20*N)))
# plot!(I2, Y2)
# plot!(p, show=true)
