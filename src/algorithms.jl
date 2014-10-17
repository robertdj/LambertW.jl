function lambertw(x::Real, k::Int=0, prec=10e-10)
	# Check branch
	if k != 0 && k != -1
		error("The branch k must be either 0 or -1")
	end

	# Check argument
	if x < -exp(-1)
		error("Lambert's W function is only defined for input greater than -exp(-1)")
	end

	if k == -1 && x >= 0
		error("The argument of branch -1 of Lambert's W function must be negative")
	end

	# First approximation
	W = lambertwApprox(x, k)

	# Compute residual
	r = abs( x - W*exp(W) )

	# Apply Halley's method to increase precision
	while r > prec
		Wnew = W - (W*exp(W) - x) / (exp(W)*(W+1) - (W+2)*(W*exp(W)-x)/(2*W+2))

		r = abs( x - W*exp(W) )
		W = Wnew
	end

	return W
end


function lambertwApprox(x::Real, k::Int)
	if x < -exp(-1)
		error("Lambert's W function is only defined for input greater than -exp(-1)")
	end

	if k == -1
		if x >= 0
			error("Input out of range for Lambert's W function")
		end

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

function lambertw(x::Array{Float64}, k::Int=0, prec=10e-10)
	W = zeros( size(x) )

	N = length(x)

	for n = 1:N
		# ----------------------------------------------------
		# If x is out of bounds, the answer is NaN
		if x[n] < -exp(-1)
			W[n] = NaN
			continue
		end
		
		if k == -1 && x[n] >= 0
			W[n] = NaN
			continue
		end
		# ----------------------------------------------------

		W[n] = lambertw( x[n], k, prec )
	end

	return W
end

