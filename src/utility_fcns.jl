# Composition of a row matrix and a quadratic polynomial
@inline function ∘(p::QuadraticPolynomial, χ::SMatrix{1,2})
    return QuadraticForm{Float64}(p.a*χ'*χ, p.b*χ, p.c)
end


# Forms lower triangular Toeplitz matrix,
# there were no currently working implementation in existing packages
function form_upper_triangular_toeplitz_matrix(v::AbstractVector{T}) where T
	N = length(v)
	G = Matrix{eltype(v)}(undef, N, N)
	for i=1:N
		for j=i:N
			G[i,j] = v[j-i+1]
		end
	end
	return G
end
