# ----------------------------------------------------------------------------------------------------------------------
# Legendre Polynomials 
# ----------------------------------------------------------------------------------------------------------------------

"""
    legendre_coefficients(n)

Calculates the coefficients of the nth Legendre polynomial using the power
series representation of Legendre polynomials,

``P_n(x) = \\frac{1}{2^n} \\sum_{k=0}^{\\lfloor n/2 \\rfloor} (-1)^k \\binom{n}{k} \\binom{2n - 2k}{n} x^{n - 2k}``.

Coefficients are ordered from lowest to highest degree for compatibility
with `evalpoly`.
"""
function legendre_coefficients(n::Int)
    nonzero = (1/2^n) .* [(-1)^k * binomial(n, Int(k)) * binomial(Int(2n - 2k), n) for k ∈ floor(n/2):-1:0]
    coefficients = zeros(n + 1)
    if n % 2 == 0
        for i ∈ 1:length(nonzero)
            coefficients[2i-1] = nonzero[i]
        end
    else
        for i ∈ 1:length(nonzero)
            coefficients[2i] = nonzero[i]
        end
    end
    return SVector{n+1}(coefficients)
end

"""
    Legendre(n)

Representation of the nth Legendre polynomial. Available methods include:
- `n`: the polynomial degree 
- `coefficients`: the polynomial coefficients 
- `polynomial`: the polynomial function 
- `weight`: the associated Laguerre polynomial weight function, ``w(x) = 1``
- `interval`: the associated Laguerre polynomial orthogonality interval, ``[-1, 1]``
"""
struct Legendre <: AbstractLegendre 
    n::Int 
    coefficients
    polynomial
    weight
    interval

    function Legendre(n)
        coefficients = legendre_coefficients(n)
        polynomial = x -> evalpoly(x, coefficients)
        weight = x -> 1
        interval = (-1, 1)
        new(n, coefficients, polynomial, weight, interval)
    end
end

"""
    LegendreP(n, x)
    LegendreP(x; n)

Evaluates the nth Legendre polynomial at x ∈ ℝ.
"""
function LegendreP(n::Int, x::T) where {T<:Real}
    return Legendre(n).polynomial(x)
end
LegendreP(x; n) = LegendreP(n, x)

"""
    P(n, x)
    P(x; n)

Evaluates the nth Legendre polynomial at x ∈ ℝ.
"""
P(n, x) = LegendreP(n, x)
P(x; n) = LegendreP(x; n)

"""
    coefficients(P)

Determines the coefficients of a given Laguerre polynomial.
"""
coefficients(P::Legendre) = P.coefficients 

"""
    weight(P)

Determines the weight function of a given Laguerre polynomial.
"""
weight(P::Legendre) = P.weight 

"""
    interval(P)

Determines the interval of orthogonality of a given Laguerre polynomial.
"""
interval(P::Legendre) = P.interval

"""
    innerproduct(P1, P2)

Calculates the inner product of the Legendre polynomials `P1` and `P2`
using the orthogonality relation

``\\int_{-1}^{1} P_m(x) P_n(x) \\: \\mathrm{d}x = \\frac{2}{2n + 1}\\delta_{nm}``
"""
function innerproduct(P1::Legendre, P2::Legendre)
    n, m = P1.n, P2.n 
    if n != m
        return 0
    else
        return 2/(2n + 1)
    end
end

# ----------------------------------------------------------------------------------------------------------------------
# Exports
# ----------------------------------------------------------------------------------------------------------------------

export Legendre
export LegendreP
export P
export coefficients 
export weight 
export interval 
export innerproduct

# ----------------------------------------------------------------------------------------------------------------------