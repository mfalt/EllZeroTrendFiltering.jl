module DynamicApproximations

export QuadraticPolynomial, PiecewiseQuadratic, QuadraticForm
export fit_pwl_constrained, fit_pwl_reguralized
export pwq_dp_regularized, pwq_dp_constrained
#Should we export the following?
export compute_transition_costs, compute_discrete_transition_costs
export recover_optimal_index_set, recover_solution, brute_force_search
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

include(joinpath("types","QuadraticPolynomial.jl"))
include(joinpath("types","PiecewiseQuadratic.jl"))
include(joinpath("types","QuadraticForm.jl"))
include(joinpath("types","transition_costs.jl"))

include("brute_force_search.jl")


AbstractTransitionCost{T} = Union{Array{QuadraticForm{T},2}, TransitionCostDiscrete{T}, TransitionCostContinuous{T}}

include("solve.jl")
include("transition_cost_computation.jl")

snp500_data() = readdlm(joinpath(Pkg.dir("DynamicApproximations"),"examples","data","snp500.txt"))

end # module
