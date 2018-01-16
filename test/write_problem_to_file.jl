function write_problem_to_file(problem_name, g, ζ_vec, I_sols, f_sols)


outfile = joinpath(Pkg.dir("DynamicApproximations"),"test",
                    "problems",  problem_name * ".jl")
f = open(outfile, "w")

println(f, "# Automatically generated test file")
println(f, "function ", problem_name, "()\n")

print(f, "g = ")
if typeof(g) <: Vector
    println(f, "[\n", join([@sprintf "%0.15f" x for x in g], ",\n"), "]\n")
else
    println(f, g)
end

println(f, "ζ_vec = ", ζ_vec, "\n")

println(f, "I_sols = [\n", join([string(x) for x in I_sols], ",\n"), "]\n")



println(f, "f_sols = [\n", join([@sprintf "%0.15f" x for x in f_sols], ",\n"), "]\n")

println(f, "return g, ζ_vec, I_sols, f_sols\n")
println(f, "end")

close(f)

end
