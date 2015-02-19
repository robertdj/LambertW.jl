# Lambert's W function for numbers
function lambertw(x::Real, k::Int=0, prec=eps())
	# Check branch
	if k != 0 && k != -1
		error("The branch k must be either 0 or -1")
	end

	# Check argument
	if x < -exp(-1)
		return NaN
	end

	if k == -1 && x >= 0
		return NaN
	end

	# First approximation
	W = lambertwApprox(x, Val{k})

	# Compute residual using logarithms to avoid numerical overflow
	# When x == 0, r = NaN, but here the approximation is exact and 
	# the while loop below is unnecessary
	r = abs( W - log(abs(x)) + log(abs(W)) )

	# Apply Fritsch’s method to increase precision
	n = 0
	while r > prec && n < 5
		n += 1

		z = log(x/W) - W
		q = 2*(1 + W)*(1 + W + 2/3*z)
		epsilon = z*(q - z) / ((1 + W)*(q - 2*z))
		W *= 1 + epsilon

		r = abs( W - log(abs(x)) + log(abs(W)) )
	end

	return W
end

# Lambert's W function for arrays
function lambertw{T<:Real}(x::Array{T}, k::Int=0, prec=eps())
	W = map( x -> lambertw(x, k, prec), x )
end


# Initial approximation for branch 0
function lambertwApprox(x::Real, ::Type{Val{0}})
	if x <= 1
		eta = 2 + 2*exp(1)*x;
		N2 = 3*sqrt(2) + 6 - (((2237+1457*sqrt(2))*exp(1) - 4108*sqrt(2) - 5764)*sqrt(eta))/((215+199*sqrt(2))*exp(1) - 430*sqrt(2)-796);
		N1 = (1-1/sqrt(2))*(N2+sqrt(2));

		W = -1 + sqrt(eta)/(1 + N1*sqrt(eta)/(N2 + sqrt(eta)));
	else
		W = log( 6*x/(5*log( 12/5*(x/log(1+12*x/5)) )) )
	end

	return W
end

# Initial approximation for branch -1
function lambertwApprox(x::Real, ::Type{Val{-1}})
	const M1 = 0.3361
	const M2 = -0.0042
	const M3 = -0.0201

	sigma = -1 - log(-x)
	W = -1 - sigma - 2/M1*(1 - 1/(1 + (M1*sqrt(sigma/2)) / (1 + M2*sigma*exp(M3*sqrt(sigma)))))

	return W
end

