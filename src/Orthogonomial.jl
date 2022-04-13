module Orthogonomial

# Imports
using Roots: find_zeros
using StaticArrays: SVector
using SpecialFunctions: gamma 

# Abstract types 
abstract type OrthogonalPolynomial end 
abstract type AbstractHermite <: OrthogonalPolynomial end 
abstract type AbstractLaguerre <: OrthogonalPolynomial end 
abstract type AbstractLegendre <: OrthogonalPolynomial end 
abstract type AbstractChebyshev <: OrthogonalPolynomial end

# Files
include("chebyshev.jl")
include("hermite.jl")
include("laguerre.jl")
include("legendre.jl")
include("utils.jl")

end