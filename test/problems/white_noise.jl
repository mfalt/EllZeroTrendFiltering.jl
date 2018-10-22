# Automatically generated test file
function white_noise()

g = [
0.297287984535462,
0.382395967790608,
-0.597634476728231,
-0.010445244637376,
-0.839026854388764,
0.311111338498334,
2.295087823837310,
-2.267086348800531,
0.529965576166746,
0.431421526422912,
0.583708287568779,
0.963271605038191,
0.458790955053717,
-0.522336757421508,
0.408395838324752]

ζ_vec = [13.2605, 5.72288, 2.46984, 1.06592, 0.46002, 0.198532, 0.0856812, 0.0369777, 0.0159586, 0.00688728]

I_sols = [
[1, 15],
[1, 3, 15],
[1, 7, 8, 15],
[1, 7, 8, 9, 15],
[1, 5, 7, 8, 9, 15],
[1, 5, 7, 8, 9, 12, 15],
[1, 5, 7, 8, 9, 12, 14, 15],
[1, 5, 6, 7, 8, 9, 12, 14, 15],
[1, 3, 4, 5, 7, 8, 9, 12, 14, 15],
[1, 2, 3, 4, 5, 7, 8, 9, 12, 14, 15],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 14, 15],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 15],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]]

f_sols = [
13.604680464796299,
13.260506698294682,
11.380933016032856,
6.070779395384360,
1.658799820085756,
1.459921430168823,
0.760578522235333,
0.586871209768404,
0.424250952535800,
0.235164307482329,
0.119283257822485,
0.048351419645437,
0.008609105509041,
-0.000000000000007]

return g, ζ_vec, I_sols, f_sols

end