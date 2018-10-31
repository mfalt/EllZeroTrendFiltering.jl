# PiecewiseLinearContinuousFitting
# Piecewise Linear Continouous Interpolation Constraints
# Sparse L2 Optimal Fitting subject Continuity Constraints
# Integrated Square

using Test

# TODO what is this?
#@everywhere include(joinpath(dirname(@__FILE__),"..","src","jl"))


@everywhere using Polynomials, IterTools#, Plots

testN = 9
@everywhere t1 = time()
#maxlengths = Array{Int64,1}(undef, testN) ;        maxlengths[1] = 1
#times = Array{Float64,1}(undef, testN);            times[1] = 0.
#sizes = Array{Array{Int64,2},1}(undef, testN);     sizes[1] = reshape([1],(1,1))
#gr()
j = 9
#@everywhere function testi(j)
    N = j
    t = range(0, stop=2π, length=50)
    #g = Poly([0,0,1,-1])
    ran = 0.1*randn(N)
    z(t) = ran[floor(Int64, t/2pi*(N-1))+1]
    g = sin
    g2(t) = g(t) + z(t)
    l = compute_transition_costs(g2, t);
    K = 10
    #@time I, Y, f = brute_force_search(l, K-1);
    Λ_0 = [create_new_pwq(minimize_wrt_x2(l[i, N])) for i in 1:N-1];
    Λ, t, _, _, _ = @timed pwq_dp_constrained(Λ_0, l, N-1);
    println("$j, tot: $(time()-t1), $t")
    times = t
    sizes = reshape([isassigned(Λ, i) ? length(Λ[i]) : 0 for i in eachindex(Λ)], size(Λ))
    maxlengths = maximum(sizes)
    #return times, sizes, maxlengths
#end
@time result = pmap(testi, 300:testN)
times, sizes, maxlengths = collect(zip(result...))
times = [times...]; sizes=[sizes...]; maxlengths=[maxlengths...];

function plotall(Λ)
    M,N = size(Λ)
    p = plot(layout=(M-2,N-1), size=(1400,1000))
    k = 1
    for m = (M-2):-1:1
        for n = 1:N-1
            idx = (n-1)*M+m
            if isassigned(Λ,idx)
                x, y, x_all, y_all = get_vals(Λ[m,n])
                for i in eachindex(x_all)
                    plot!(p, y_all[i][[1,end]], x_all[i][[1,end]], l = 0, m=(1,:cross,:orange), subplot=k, lab="")
                    plot!(p, y_all[i], x_all[i], l=(1,:blue), lab="", subplot=k)
                end
                plot!(p, y, x, l=:red, subplot=k, lab="$m, $n, $k", ylims=(-4,4), xlims=(0,6))
            end
            k += 1
        end
    end
    return p
end

plot(times)
plot!(collect(500:512),times3*2)
plot!(((1:512).^3.3)*maximum(times)/300^3.3)
plot!(ylims=(0,1400))
plot!()
heatmap!(sizes[end])

plot!(maximum(times)*sum.(sizes)/maximum(sum.(sizes)))

using Plots
plot(t, g.(t), lab="g(t)")
for k = 7:7
    @time I2, y2, f2 = recover_solution(Λ[k, 1], 1, N)
    Y2, _ = find_optimal_y_values(l, I2)

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
