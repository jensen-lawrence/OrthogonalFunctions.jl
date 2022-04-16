# Functions
"""
    degree(p)

Returns the degree of `p`.
"""
degree(p::AbstractOrthogonalPolynomial) = p.n 

"""
    polynomial(p)

Returns the polynomial function of `p`.
"""
polynomial(p::AbstractOrthogonalPolynomial) = p.f 

# Files 
include("hermite/hermite.jl")
# include("jacobi/jacobi.jl")
# include("laguerre/laguerre.jl")

# Exports 
export degree
export polynomial