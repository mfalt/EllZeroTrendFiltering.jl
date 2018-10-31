#### This example shows the growth of the number of quadratic functions in the constrained case
#### The figure in the article was generated with "OPTIMIZE = false"

using EllZeroTrendFiltering
using DSP
using Plots
using Interpolations

### START TEST FUNCTIONS


t = 1:200
srand(1)

circle_segment(N) = sin.(acos.(range(-1, stop=1, length=N)))
linear_trend(N) = (0:N-1)

g_cll = [0.2*linear_trend(120) + circle_segment(120);
     -6 + 8*linear_trend(40);
     -6 + 8*linear_trend(40)]

g_cllc = [5*circle_segment(70);
    -5 + 4*linear_trend(30);
    -5 + 4*linear_trend(30);
     5*circle_segment(70)]


g_rand = () -> randn(200)

g_walk = () -> cumsum(randn(200))

g_square = range(0, stop=4, length=200).^2

g_cube = range(0, stop=4, length=200).^3

g_4 = range(0, stop=4, length=200).^4

g_5 = range(0, stop=4, length=200).^5

g_exp = exp.(range(0, stop=4, length=200))

g_super_exp = exp.(range(0, stop=3, length=200).^2)

g_super_exp2 = [exp.(range(0, stop=2, length=100).^2)   ;  exp.(range(2, stop=0.5, length=100).^2)]

g_super_exp3 = [exp.(range(0, stop=2, length=100).^2)   ; -exp.(range(2, stop=0.5, length=100).^2)]

g_super_exp4 = [exp.(range(2, stop=0.5, length=100).^2) ;  exp.(range(0, stop=2, length=100).^2)]

g_super_exp5 = [exp.(range(2, stop=0, length=100).^2)   ; 2-exp.(range(0, stop=0, length=100).^2)]

function randlin(N)
    t_grid = cumsum( rand(10:25, N) )
    y_grid = randn(N)
    g = interpolate((t_grid,), y_grid, Gridded(Linear()))[t]
    reverse!(g)
end

function randlinfilt1(N)
    t_grid = cumsum( rand(10:25, N) )
    y_grid = randn(N)
    g = interpolate((t_grid,), y_grid, Gridded(Linear()))[t]
    g = reverse(filt(ones(10), 1, reverse(filt(ones(10), 1, g))))
    reverse!(g)
end

function randlinfilt2(N)
    t_grid = cumsum( rand(10:25, N) )
    y_grid = randn(N)
    g = interpolate((t_grid,), y_grid, Gridded(Linear()))[t]
    g =filt(ones(40), 1, g)
end

function randlinfilt3(N)
    t_grid = cumsum( rand(10:25, N) )
    y_grid = randn(N)
    g = interpolate((t_grid,), y_grid, Gridded(Linear()))[t]
    g = filtfilt(ones(10), g)
end


g_data = [g_cll, g_cllc,
        g_exp, g_super_exp, g_super_exp2,
        g_super_exp3, g_super_exp4, g_super_exp5,
        g_square, g_cube, g_4, g_5,
        reverse!(g_cll), reverse!(g_cllc),
        reverse!(g_exp), reverse!(g_super_exp), reverse!(g_super_exp2),
        reverse!(g_super_exp3), reverse!(g_super_exp4), reverse!(g_super_exp5),
        reverse!(g_square), reverse!(g_cube), reverse!(g_4), reverse!(g_5),
        g_rand(), g_rand(), g_rand(), g_rand(), g_rand(), g_rand(), g_rand(), g_rand(), g_rand(), g_rand(),
        g_walk(), g_walk(), g_walk(), g_walk(), g_walk(), g_walk(), g_walk(), g_walk(), g_walk(), g_walk(),
        randlin(5), randlinfilt1(5), randlinfilt2(5), randlinfilt3(5),
        randlin(10), randlinfilt1(10), randlinfilt2(10), randlinfilt3(10),
        randlin(20), randlinfilt1(20), randlinfilt2(20), randlinfilt3(20),
        randlin(40), randlinfilt1(40), randlinfilt2(40), randlinfilt3(40),
        reverse!(randlin(5)), reverse!(randlinfilt1(5)), reverse!(randlinfilt2(5)), reverse!(randlinfilt3(5)),
        reverse!(randlin(10)), reverse!(randlinfilt1(10)), reverse!(randlinfilt2(10)), reverse!(randlinfilt3(10)),
        reverse!(randlin(20)), reverse!(randlinfilt1(20)), reverse!(randlinfilt2(20)), reverse!(randlinfilt3(20)),
        reverse!(randlin(40)), reverse!(randlinfilt1(40)), reverse!(randlinfilt2(40)), reverse!(randlinfilt3(40)),
        cumsum(randlin(5)), cumsum(randlinfilt1(5)), cumsum(randlinfilt2(5)), cumsum(randlinfilt3(5)),
        cumsum(randlin(10)), cumsum(randlinfilt1(10)), cumsum(randlinfilt2(10)), cumsum(randlinfilt3(10)),
        cumsum(randlin(20)), cumsum(randlinfilt1(20)), cumsum(randlinfilt2(20)), cumsum(randlinfilt3(20)),
        cumsum(randlin(40)), cumsum(randlinfilt1(40)), cumsum(randlinfilt2(40)), cumsum(randlinfilt3(40)),
        ]
### END TEST FUNCTIONS

nbr_prob = length(g_data)

V = zeros(length(t), nbr_prob)

g_save = zeros(length(t), nbr_prob)

for prob_i=1:nbr_prob
    g = g_data[prob_i]

    l = compute_discrete_transition_costs(g, t);
    V_N = QuadraticPolynomial(0.0, 0.0, 0.0)

    Λ =  pwq_dp_constrained(l, V_N, 15);

    l = [(isassigned(Λ,m,i) ? length(Λ[m,i]) : 0) for m = 1:size(Λ,1), i = 1:size(Λ,2)]

    V[:, prob_i] .= maximum(l, 1)[:]
    g_save[:, prob_i] .= g
end



gr()
plot(size=(800,400));
plot(t, V, lab="", linealpha=0.2, c=:black);
plot!(t, length(t) - t + 1, l=(:red, :dash), lab="")
savefig("complexity.png")
writecsv("complexity.csv", V)
