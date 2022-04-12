module Orthogonomial

# Imports
using StaticArrays: SVector
using SpecialFunctions: gamma 

# Abstract types 
abstract type OrthogonalPolynomial end 
abstract type AbstractHermite <: OrthogonalPolynomial end 
abstract type AbstractLaguerre <: OrthogonalPolynomial end 
abstract type AbstractJacobi <: OrthogonalPolynomial end 
abstract type AbstractChebyshev <: AbstractJacobi end 
abstract type AbstractLegendre <: AbstractJacobi end

# Files
include("chebyshev.jl")
include("hermite.jl")
include("laguerre.jl")
include("legendre.jl")
include("utils.jl")

end