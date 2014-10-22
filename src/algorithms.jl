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

	# Compute function
	W = lambertwNAC(x, k, prec)

	return W
end


# NAC: No Argument Check
function lambertwNAC(x::Real, k::Int=0, prec=eps())
	# First approximation
	W = lambertwApprox(x, k)

	# Compute residual using logarithms to avoid numerical overflow
	r = abs( W - log(abs(x)) + log(abs(W)) )

	# Apply Halley's method to increase precision
	n = 0
	while r > prec && n < 5
		n += 1

		Wnew = W - (W*exp(W) - x) / (exp(W)*(W+1) - (W+2)*(W*exp(W)-x)/(2*W+2))

		r = abs( W - log(abs(x)) + log(abs(W)) )
		W = Wnew
	end

	return W
end


function lambertwApprox(x::Real, k::Int)
	if k == -1
		M1 = 0.3361
		M2 = -0.0042
		M3 = -0.0201

		sigma = -1 - log(-x)
        W = -1 - sigma - 2/M1*(1 - 1/(1 + (M1*sqrt(sigma/2)) / (1 + M2*sigma*exp(M3*sqrt(sigma)))))
	elseif k == 0
		if x <= 1
            eta = 2 + 2*exp(1)*x;
            N2 = 3*sqrt(2) + 6 - (((2237+1457*sqrt(2))*exp(1) - 4108*sqrt(2) - 5764)*sqrt(eta))/((215+199*sqrt(2))*exp(1) - 430*sqrt(2)-796);
            N1 = (1-1/sqrt(2))*(N2+sqrt(2));

            W = -1 + sqrt(eta)/(1 + N1*sqrt(eta)/(N2 + sqrt(eta)));
		else
            W = log( 6*x/(5*log( 12/5*(x/log(1+12*x/5)) )) )
		end
	end

	return W
end


# ------------------------------------------------------------
# Lambert's W function for arrays
# ------------------------------------------------------------

function lambertw{T<:Real}(x::Array{T}, k::Int=0, prec=eps())
	# Check branch
	if k != 0 && k != -1
		error("The branch k must be either 0 or -1")
	end

	# Constants
	lower_bound = -exp(-1)
	lower_branch = k == -1

	W = zeros( size(x) )

	N = length(x)

	for n = 1:N
		# ----------------------------------------------------
		# If x is out of bounds, the answer is NaN
		if x[n] < lower_bound
			W[n] = NaN
			continue
		end
		
		if lower_branch && x[n] >= 0
			W[n] = NaN
			continue
		end
		# ----------------------------------------------------

		W[n] = lambertwNAC( x[n], k, prec )
	end

	return W
end

