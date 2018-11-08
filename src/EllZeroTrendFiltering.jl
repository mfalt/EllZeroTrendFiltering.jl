module EllZeroTrendFiltering

export QuadraticPolynomial, PiecewiseQuadratic, QuadraticForm
export fit_pwl_constrained, fit_pwl_regularized
export construct_value_fcn_constrained, construct_value_fcn_regularized
#Should we export the following?
export compute_problem_data_pwl, compute_problem_data_lti
export recover_optimal_index_set, recover_solution, brute_force_search
export recover_optimal_index_set_zero_ic
export generate_markov_matrix
export snp500_data

import Base.-
import Base.+
import Base.∘
import Base.show
import Base: iterate, length, zero, getindex, ==

import IterTools
import ControlSystems
using StaticArrays
using QuadGK
using MatrixDepot

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


include("utility_fcns.jl")
include("insert_quadratic_poly.jl")
include("compute_problem_data.jl")
include("brute_force_search.jl")
include("solve.jl")


snp500_data() = readdlm(joinpath(dirname(@__FILE__),"..","examples","data","snp500.txt"))

end # module
