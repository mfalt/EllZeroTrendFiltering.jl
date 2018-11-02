"""
    l = precompute_transition_costs(l::AbstractTransitionCost{T})
Precompute the transition costs `l` and returns a matrix of quadratic forms.
This might be useful if performing multiple computations with the same problem data.
"""
function precompute_transition_costs(l::AbstractTransitionCost{T}) where T
    N = size(l, 2)
    l_precomputed = Array{QuadraticForm{T}}(undef, N-1,N)
    for i = 1:N-1
        for ip = i+1:N
            l_precomputed[i,ip] = l[i,ip]
        end
    end
    return l_precomputed
end


# Continuous case
""" `compute_problem_data_pwl(g, t, lazy; tol=1e-3)`
    Return `(l, V_N)`: the transition cost and end-cost.
    Slightly unsafe since `tol` is ignored (not needed) on discrete problems
"""
function compute_problem_data_pwl(g, t; precompute=false, tol=1e-3)
    T = Float64

    l = TransitionCostContinuous{T}(g, t, tol)
    if precompute
        l = precompute_transition_costs(l)
    end

    χ = form_upper_triangular_toeplitz_matrix([(@SMatrix [0 one(T)]) for k=1:length(t)])

    # Continouous case, no cost at endpoint
    V_N = zero(QuadraticPolynomial{Float64})
    return l, χ, V_N
end

# Discrete case, tol is not used here, but in signature to enable dispatch
function compute_problem_data_pwl(g::AbstractArray, t=1:length(g); precompute=false, tol=1e-3)
    T = promote_type(eltype(g), Float64)

    l = TransitionCostDiscrete{T}(g, t)
    if precompute
        l = precompute_transition_costs(l)
    end

    χ = form_upper_triangular_toeplitz_matrix([(@SMatrix [0 one(T)]) for k=1:length(t)])

    # Discrete case, so cost at endpoint is quadratic
    V_N = QuadraticPolynomial(1.0, -2*g[t[end]], g[t[end]]^2)
    return l, χ, V_N
end



# As of now, it is implicitly understood that B = [0; 1]
function compute_problem_data_lti(g::AbstractArray, A::AbstractMatrix, C::AbstractMatrix)

	T = Float64
	A = SMatrix{2,2,T,4}(A)
	C = SMatrix{1,2,T,2}(C)

	N = length(g)

	Apow = Vector{SMatrix{2,2,T,4}}(undef, N)
	for k=1:N
		Apow[k] = A^(k-1)
	end

	l = Matrix{QuadraticForm{T}}(undef, N, N)
	for i=1:N-1
		l[i,i] = QuadraticForm(SMatrix{2,2,T,4}(zeros(4)), SVector{2,T}(zeros(2)), zero(T))

	    for j=i:N-1
	    	# Need parenthesis on A to avoid chnage of type tp just Array
	        l[i,j+1] = QuadraticForm(
				        	l[i,j].P + Apow[j-i+1]'*C'*C*Apow[j-i+1],
				        	l[i,j].q - (2g[j]*C*Apow[j-i+1])[:], # FIXME problem with getting matrix/vector right
				        	l[i,j].r + g[j]^2)
	    end
	end

	χ = form_upper_triangular_toeplitz_matrix([(@SMatrix [1.0 0]) * Apow[k] for k=1:N])

	V_N = QuadraticPolynomial(1.0, -2*g[end], g[end]^2)

	return l, χ, V_N
end
