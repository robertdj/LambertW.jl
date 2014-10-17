using LambertW
using Base.Test

# Test special values
@test lambertw(-exp(-1), 0) == -1
@test lambertw(-exp(-1), -1) == -1

@test lambertw(0) == 0

a = linspace( exp(-1), exp(1), 10 )
for n = 1:10
	@test_approx_eq lambertw( -log(a[n])/a[n] ) -log(a[n])
end

# The r.h.s. is the Omega constant
@test_approx_eq lambertw(1) 0.5671432904097838729999686622

