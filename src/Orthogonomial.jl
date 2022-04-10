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
include("hermite.jl")

end