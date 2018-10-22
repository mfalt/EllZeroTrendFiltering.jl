# Automatically generated test file
function super_exponential3()

g = [
1.000000000000000,
1.064494458917860,
1.284025416687741,
1.755054656960298,
2.718281828459045,
4.770733181967603,
9.487735836358526,
21.380942759123343,
54.598150033144236,
-54.598150033144236,
-21.380942759123343,
-9.487735836358526,
-4.770733181967603,
-2.718281828459045,
-1.755054656960298,
-1.284025416687741]

ζ_vec = [6649.18, 1320.94, 262.421, 52.1332, 10.3569, 2.05753, 0.408753, 0.0812037, 0.0161321, 0.00320484]

I_sols = [
[1, 16],
[1, 8, 16],
[1, 9, 10, 16],
[1, 7, 9, 10, 16],
[1, 7, 9, 10, 11, 16],
[1, 6, 8, 9, 10, 11, 16],
[1, 6, 8, 9, 10, 11, 12, 16],
[1, 5, 7, 8, 9, 10, 11, 12, 16],
[1, 5, 7, 8, 9, 10, 11, 12, 13, 16],
[1, 4, 6, 7, 8, 9, 10, 11, 12, 13, 16],
[1, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16],
[1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16],
[1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
[1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]]

f_sols = [
6906.314103950247045,
6649.182023111979106,
1780.147768604715566,
810.634794153323128,
149.682164816596014,
79.779903677661423,
16.245588913173378,
8.328684761258955,
1.929086519737211,
0.854273036134146,
0.251603511120265,
0.081666565196429,
0.041290097946330,
0.004006052640761,
-0.000000000020918]

return g, ζ_vec, I_sols, f_sols

end