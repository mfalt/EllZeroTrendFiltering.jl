# Automatically generated test file
function discontinuous1()

g = [
0.000000000000000,
0.828539361054709,
1.231479419283098,
1.542809041582063,
1.793807989999907,
1.993807989999906,
2.142809041582063,
2.231479419283098,
2.228539361054709,
1.800000000000000,
-6.000000000000000,
2.000000000000000,
10.000000000000000,
18.000000000000000,
26.000000000000000,
-6.000000000000000,
2.000000000000000,
10.000000000000000,
18.000000000000000,
26.000000000000000]

ζ_vec = [873.808, 250.088, 71.5766, 20.4856, 5.86308, 1.67804, 0.480265, 0.137454, 0.0393402, 0.0112594]

I_sols = [
[1, 20],
[1, 18, 20],
[1, 15, 16, 20],
[1, 12, 15, 16, 20],
[1, 10, 11, 15, 16, 20],
[1, 5, 10, 11, 15, 16, 20],
[1, 3, 8, 10, 11, 15, 16, 20],
Int[],
[1, 2, 4, 7, 9, 10, 11, 15, 16, 20]]

f_sols = [
1102.237767513050130,
873.808100414683622,
544.329717092552983,
57.308828172005178,
1.311071990397977,
0.282505163300357,
0.069492548637299,
0.018264420672040,
0.004190223377918]

return g, ζ_vec, I_sols, f_sols

end