# LambertW #

[![Build Status](https://travis-ci.org/git@dahl-jacobsen.dk/LambertW.jl.svg?branch=master)](https://travis-ci.org/git@dahl-jacobsen.dk/LambertW.jl)

This package implements an approximation of the two real branches of [Lambert's W function](https://en.wikipedia.org/wiki/Lambert_W_function), W = W(x), defined as a solution of the equation W*exp(W) = x

The approximation is based on an iterative method.


## Usage ##

To compute W_k(x), where k being either 0 or -1 denotes the branch, use
	
	lambertw(x, k)

A third, optional argument specifies the precision.


## References ##

The initial approximations are described in

* D.A. Barry, J.-Y. Parlange, L. Li, H. Prommer, C. J. Cunningham, F. Stagnitti, "Analytical approximations for real values of the Lambert W-function", Mathematics and Computers in Simulation, volume 53, 2000, page 95--103.
DOI: [10.1016/S0378-4754(00)00172-5](https://dx.doi.org/10.1016/S0378-4754(00)00172-5)
* The erratum to the above article; DOI: [10.1016/S0378-4754(02)00051-4](https://dx.doi.org/10.1016/S0378-4754(02)00051-4)

The iterative refinement is calculated using Halley's method as mentioned in e.g.

* Darko Veberic, "Lambert W function for applications in physics", Computer Physics Communications, volume 183, 2012, page 2622-–2628.
DOI: [10.1016/j.cpc.2012.07.008](http://dx.doi.org/10.1016/j.cpc.2012.07.008)

