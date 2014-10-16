function lambertw(x::Real, k::Int=0, prec=10e-5)
	# Check input
	if k != 0 && k != -1
		error("The branch k must be either 0 or -1")
	end

	if x < -exp(-1)
		error("Lambert's W function is only defined for input greater than -exp(-1)")
	end

	if k == -1 && x >= 0
		error("The argument of branch -1 of Lambert's W function must be negative")
	end

	# TODO: Include special values

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
	# TODO: Find better approx
	if x > 1
		logx = log(x)
		llogx = log(log(x))

		W = logx - llogx + llogx/logx
	else
		W = 0
	end

	return W
end

