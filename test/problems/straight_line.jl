function straight_line()

    g = 0:0.02:1

    M_bf = 2

    ζ_vec = [100, 30, 10, 3, 1, 0.3, 0.1, 0.03, 0.01]

    # For each m there are many non-unique solution
    I_sols =
    [[[1, length(g)]];
     [Vector{Int}() for m=1:9]]

    f_sols = zeros(Float64, 10)

    return g, ζ_vec, I_sols, f_sols
end
