using Test
import EllZeroTrendFiltering: TransitionCostContinuous

#TODO let to avoid global scope warning, change back in julia 1.0?
#let
t = [-1; 0; 1]

g = x -> 0
l = TransitionCostContinuous{Float64}(g, t)
@test l[1,2].P == 1/6*[2 1; 1 2]
@test l[1,2].q == [0,0]
@test l[1,2].r == 0


g = x -> 1
l = TransitionCostContinuous{Float64}(g, t)
@test l[1,2].P == 1/6*[2 1; 1 2]
@test l[1,2].q == [-1,-1]
@test l[1,2].r == 1


g = x -> x
l = TransitionCostContinuous{Float64}(g, t)
@test l[2,3].P == 1/6*[2 1; 1 2]
@test l[2,3].q ≈ [-1/3,-2/3]
@test l[2,3].r ≈ 1/3
@test l[2,3](0,0) == 1/3
@test l[2,3](1,1) == 1/3
@test l[2,3](1,2) ≈ 1
@test l[2,3](1,0) ≈ 1/3
#end
