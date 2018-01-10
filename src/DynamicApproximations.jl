module DynamicApproximations

export QuadraticPolynomial, PiecewiseQuadratic, QuadraticForm
export fit_pwl_constrained, fit_pwl_reguralized
export regularize, find_optimal_fit
#Should we export the following?
export compute_transition_costs, recover_solution, brute_force_optimization
export recover_optimal_index_set, compute_discrete_transition_costs
export snp500_data

import Base.-
import Base.+
import Base.show
import Base.Operators
import Base: start, next, done, length, zero, getindex, ==

import IterTools

using StaticArrays
using QuadGK

global const DEBUG = false
global const DEBUG2 = false
global const COUNTER_TEST = false
global const OPTIMIZE = true
# TODO Tweak this accuracy, sqrt(eps()) was not good enough for snp 500, M = 12 test
global const accuracy = sqrt(eps())/100000


include(joinpath("types","QuadraticPolynomial.jl"))
include(joinpath("types","PiecewiseQuadratic.jl"))
include(joinpath("types","QuadraticForm.jl"))
include(joinpath("types","TransitionCostLazy.jl"))


AbstractTransitionCost{T} = Union{Array{QuadraticForm{T},2}, AbstractTransitionCostLazy{T}}

include("solve.jl")
include("transition_cost_computation.jl")

snp500_data() = readdlm(joinpath(Pkg.dir("DynamicApproximations"),"examples","data","snp500.txt"))

end # module
