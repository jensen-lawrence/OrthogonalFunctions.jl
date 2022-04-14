module Orthogonomial

# Imports
using Roots: find_zeros
using SpecialFunctions: gamma 

# Abstract types 
abstract type OrthogonalPolynomial end 
abstract type AbstractHermite <: OrthogonalPolynomial end 
abstract type AbstractLaguerre <: OrthogonalPolynomial end 
abstract type AbstractJacobi <: OrthogonalPolynomial end
abstract type AbstractChebyshev <: AbstractJacobi end

# Files
include("chebyshev.jl")
include("gegenbauer.jl")
include("hermite.jl")
include("jacobi.jl")
include("laguerre.jl")
include("legendre.jl")
include("utils.jl")

# Exports 
export degree, polynomial, weight, interval
export innerproduct, roots 
export ChebyshevT, T, ChebyshevFirst
export ChebyshevU, U, ChebyshevSecond 
export ChebyshevV, V, ChebyshevThird 
export ChebyshevW, W, ChebyshevFourth 
export GegenbauerC, C, Gegenbauer, Ultraspherical 
export JacobiP, P, Jacobi 
export LegendreP, P, Legendre

end