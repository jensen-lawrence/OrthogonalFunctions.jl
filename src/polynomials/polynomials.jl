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

"""
    roots(p)

Returns the roots of `p`.
"""
function roots(p::AbstractOrthogonalPolynomial)
    @vars x 
    return nroots(p(x))
end

# Files 
include("hermite/hermite.jl")
# include("jacobi/jacobi.jl")
# include("laguerre/laguerre.jl")

# Exports 
export degree
export polynomial
export roots