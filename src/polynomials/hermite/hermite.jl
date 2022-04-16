# Functions
"""
    interval(H)

Returns the orthogonality interval of H.
"""
function interval(H::AbstractHermitePolynomial)
    return (-Inf, Inf)
end

"""
    roots(H)

Returns the roots (zeros) of H.
"""
function roots(H::AbstractHermitePolynomial)
    bound = 2*√(H.n)
    return find_zeros(H, (-bound, bound))
end

# Files 
include("hermite_polynomials.jl")
include("probabilistic_hermite_polynomials.jl")

# Exports 
export HermitePolynomial 
export HermiteH 
export H 
export ProbabilisticHermitePolynomial
export HermiteHe 
export He 
export weight 
export innerproduct 
export derivative 
export interval 
export roots