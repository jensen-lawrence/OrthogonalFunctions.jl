# Functions
"""
    interval(H)

Returns the orthogonality interval of H.
"""
function interval(H::AbstractHermitePolynomial)
    return (-Inf, Inf)
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