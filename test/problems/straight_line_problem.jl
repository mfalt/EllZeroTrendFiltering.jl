function straight_line_problem()

    g = 0:0.02:1

    M_bf = 2

    ζvec = [100, 30, 10, 3, 1, 0.3, 0.1, 0.03, 0.01]

    # For each m there are many non-unique solution
    I_sols =
    [[[1, length(g)]];
     [Vector{Int}() for m=1:9]]

    f_sols = zeros(Float64, 10)

    return g, M_bf, ζvec, I_sols, f_sols
end
