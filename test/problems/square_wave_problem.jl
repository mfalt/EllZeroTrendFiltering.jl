function square_wave_problem()

    g = kron([1, -1, 1, -1, 1], ones(5))

    M_bf = 5

    ζvec = [100, 30, 10, 3, 1, 0.3, 0.1, 0.03, 0.01]

    # For each m there are many non-unique solution
    I_sols =
    [Vector{Int}() for m=1:24]

    f_sols =
        [[24.0,
        18.10194908759502,
        15.757575757575758,
        6.602807597027249,
        4.84848484848485,
        4.285714285714288,
        2.4242424242424248,
        2.1356521739130434];
        zeros(Float64, 16)]

    return g, M_bf, ζvec, I_sols, f_sols
end
