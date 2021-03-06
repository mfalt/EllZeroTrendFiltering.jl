module EllZeroTrendFiltering

export QuadraticPolynomial, PiecewiseQuadratic, QuadraticForm
export fit_pwl_constrained, fit_pwl_regularized
export pwq_dp_regularized, pwq_dp_constrained
#Should we export the following?
export compute_transition_costs, compute_discrete_transition_costs
export recover_optimal_index_set, recover_solution, brute_force_search
export snp500_data

import Base.-
import Base.+
import Base.show
import Base: iterate, length, zero, getindex, ==

import IterTools



using StaticArrays
using QuadGK

#StdLibrary
using LinearAlgebra
import DelimitedFiles: readdlm
import Printf: @printf
#Pkg, DelimitedFiles
global const DEBUG = false
global const DEBUG2 = false
global const COUNTER_TEST = false
global const OPTIMIZE = true

# TODO Tweak this accuracy, sqrt(eps()) was not good enough for snp 500, M = 12 test
global const ACCURACY = sqrt(eps())/100000

include(joinpath("types","QuadraticPolynomial.jl"))
include(joinpath("types","PiecewiseQuadratic.jl"))
include(joinpath("types","QuadraticForm.jl"))
include(joinpath("types","transition_costs.jl"))

AbstractTransitionCost{T} = Union{Array{QuadraticForm{T},2}, TransitionCostDiscrete{T}, TransitionCostContinuous{T}}

include("transition_cost_computation.jl")
include("brute_force_search.jl")
include("solve.jl")


snp500_data() = readdlm(joinpath(dirname(@__FILE__),"..","examples","data","snp500.txt"))

end # module
