module LambertW

# package code goes here
export
    lambertw

end # module


# Lambert's W function for numbers
@doc """
```
lambertw(x[, k, prec])
```

The `k`'th branch of Lambert's W function evaluated at `x`.

`k` must be `-1` or `0` (default).
`prec` defaults to `eps()`.
"""->
function lambertw(x::Real, k::Int=0, prec=eps())
	# Check branch
	if k != 0 && k != -1
		error("The branch k must be either 0 or -1")
	end

	# Check argument
	const minx = -1/e
	if x < minx
		return NaN
	end

	if k == -1 && x >= 0
		return NaN
	end

	# First approximation
	W = lambertw(x, Val{k})

	# Compute residual using logarithms to avoid numerical overflow
	# When x == 0, r = NaN, but here the approximation is exact and 
	# the while loop below is unnecessary
	r = abs( W - log(abs(x)) + log(abs(W)) )

	# Apply Fritsch’s method to increase precision
	n = 1
	while r > prec && n <= 5
		z = log(x/W) - W
		q = 2*(1 + W)*(1 + W + 2/3*z)
		epsilon = z*(q - z) / ((1 + W)*(q - 2*z))
		W *= 1 + epsilon

		r = abs( W - log(abs(x)) + log(abs(W)) )
		n += 1
	end

	return W
end

# Lambert's W function for arrays
function lambertw{T<:Real}(x::Array{T}, k::Int=0, prec=eps())
	W = map( x -> lambertw(x, k, prec), x )
end


# ------------------------------------------------------------
# Initial approximations

# For branch 0
function lambertw(x::Real, ::Type{Val{0}})
	if x <= 1
		const sqrt2 = sqrt(2)

		eta = 2 + 2*e*x;
		sqeta = sqrt(eta)

		N2 = 3*sqrt2 + 6 - (((2237+1457*sqrt2)*e - 4108*sqrt2 - 5764)*sqeta)/((215+199*sqrt2)*e - 430*sqrt2-796);
		N1 = (1-1/sqrt2)*(N2+sqrt2);

		W = -1 + sqeta/(1 + N1*sqeta/(N2 + sqeta));
	else
		W = log( 6*x/(5*log( 12/5*(x/log(1+12*x/5)) )) )
	end

	return W
end

# For branch -1
function lambertw(x::Real, ::Type{Val{-1}})
	const M1 = 0.3361
	const M2 = -0.0042
	const M3 = -0.0201

	sigma = -1 - log(-x)
	W = -1 - sigma - 2/M1*(1 - 1/(1 + (M1*sqrt(sigma/2)) / (1 + M2*sigma*exp(M3*sqrt(sigma)))))

	return W
end

