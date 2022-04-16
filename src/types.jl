abstract type AbstractOrthogonalFunction end 
abstract type AbstractOrthogonalPolynomial <: AbstractOrthogonalFunction end 
abstract type AbstractHermitePolynomial <: AbstractOrthogonalPolynomial end 
abstract type AbstractLaguerrePolynomial <: AbstractOrthogonalPolynomial end 
abstract type AbstractJacobiPolynomial <: AbstractOrthogonalPolynomial end 
abstract type AbstractChebyshevPolynomial <: AbstractJacobiPolynomial end