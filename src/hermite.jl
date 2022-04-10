# ----------------------------------------------------------------------------------------------------------------------
# Hermite Polynomials 
# ----------------------------------------------------------------------------------------------------------------------

"""
    hermite_coefficients(n)

Calculates the coefficients of the nth Hermite polynomial using the power
series representation of Hermite polynomials,

``H_n(x) = n! \\sum_{k=0}^{\\lfloor n/2 \\rfloor} \\frac{(-1)^k}{k!(n - 2k)!} (2x)^{n - 2k}``.

Coefficients are ordered from lowest to highest degree for compatibility
with `evalpoly`.
"""
function hermite_coefficients(n::Int)
    if n ≤ 20
        nonzero = factorial(n) .* [((-1.0)^k * 2.0^(n - 2k))/(factorial(k) * factorial(n - 2k)) for k ∈ floor(n/2):-1:0]
    else
        nonzero = gamma(n + 1) .* [((-1.0)^k * 2.0^(n - 2k))/(gamma(k + 1) * gamma(n - 2k + 1)) for k ∈ floor(n/2):-1:0]
    end

    nonzero = map(x -> round(x), nonzero)

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
    Hermite(n)

Representation of the nth Hermite polynomial. Available methods include:
- `n`: the polynomial degree 
- `coefficients`: the polynomial coefficients 
- `polynomial`: the polynomial function 
- `weight`: the Hermite polynomial weight function, ``w(x) = e^{-x^2}``
- `interval`: the Hermite polynomial orthogonality interval, ``(-\\inf, \\inf)``
"""
struct Hermite <: AbstractHermite 
    n::Int
    coefficients
    polynomial
    weight
    interval

    function Hermite(n)
        coefficients = hermite_coefficients(n)
        polynomial = x -> evalpoly(x, coefficients)
        weight = x -> exp(-x^2)
        interval = (-Inf, Inf)
        new(n, coefficients, polynomial, weight, interval)
    end
end

# Hermite polynomials functions
"""
    HermiteH(n, x)
    HermiteH(x; n)

Evaluates the nth Hermite polynomial at x ∈ ℝ.
"""
function HermiteH(n::Int, x::T) where {T<:Real}
    return Hermite(n).polynomial(x)
end
HermiteH(x; n) = HermiteH(n, x)

"""
    H(n, x)
    H(x; n)

Evaluates the nth Hermite polynomial at x ∈ ℝ.
"""
H(n, x) = HermiteH(n, x)
H(x; n) = HermiteH(x; n)

# Functions for getting properties 
"""
    coefficients(H)

Determines the coefficients of a given Hermite polynomial.
"""
coefficients(H::Hermite) = H.coefficients 

"""
    weight(H)

Determines the weight function of a given Hermite polynomial.
"""
weight(H::Hermite) = H.weight 

"""
    interval(H)

Determines the interval of orthogonality of a given Hermite polynomial.
"""
interval(H::Hermite) = H.interval 

"""
    innerproduct(H1, H2)

Calculates the inner product of the Hermite polynomials `H1` and `H2`
using the orthogonality relation 

``\\int_{-\\inf}^{\\inf} H_m(x) H_n(x) e^{-x^2} \\: \\mathrm{d}x = \\sqrt{\\pi} 2^n n! \\delta_{nm}``.
"""
function innerproduct(H1::Hermite, H2::Hermite)
    n, m = H1.n, H2.m
    if n != m
        return 0
    else
        if n ≤ 20
            n! = factorial(n)
        else
            n! = gamma(n + 1)
        end
        return √π * 2^n * n!
    end
end

# ----------------------------------------------------------------------------------------------------------------------
# Probabilist's Hermite Polynomials 
# ----------------------------------------------------------------------------------------------------------------------

"""
    probabilist_hermite_coefficients(n)

Calculates the coefficients of the nth probabilist's Hermite polynomial
using the power series representation of probabilist's Hermite polynomials,

``He_n(x) = n! \\sum_{k=0}^{\\lfloor n/2 \\rfloor} \\frac{(-1)^k}{k!(n - 2k)!} \\frac{x^{n - 2k}}{2^k}``.

Coefficients are ordered from lowest to highest degree for compatibility
with `evalpoly`.
"""
function probabilist_hermite_coefficients(n::Int)
    if n ≤ 20
        nonzero = factorial(n) .* [((-1.0)^k)/(2.0^k * factorial(k) * factorial(n - 2k)) for k ∈ floor(n/2):-1:0]
    else
        nonzero = gamma(n + 1) .* [((-1.0)^k)/(2.0^k * gamma(k + 1) * gamma(n - 2k + 1)) for k ∈ floor(n/2):-1:0]
    end

    nonzero = map(x -> round(x), nonzero)

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
    ProbabilistHermite(n)

Representation of the nth probabilist's Hermite polynomial. Available methods
include:
- `n`: the polynomial degree 
- `coefficients`: the polynomial coefficients 
- `polynomial`: the polynomial function 
- `weight`: the probabilist's Hermite polynomial weight function, ``w(x) = e^{-(x^2)/2}``
- `interval`: the probabillist's Hermite polynomial orthogonality interval, ``(-\\inf, \\inf)``
"""
struct ProbabilistHermite <: AbstractHermite 
    n::Int
    coefficients
    polynomial
    weight
    interval

    function ProbabilistHermite(n)
        coefficients = probabilist_hermite_coefficients(n)
        polynomial = x -> evalpoly(x, coefficients)
        weight = x -> exp(-(x^2)/2)
        interval = (-Inf, Inf)
        new(n, coefficients, polynomial, weight, interval)
    end
end

"""
    HermiteHe(n, x)
    HermiteHe(x; n)

Evaluates the nth probabilist's Hermite polynomial at x ∈ ℝ.
"""
function HermiteHe(n::Int, x::T) where {T<:Real}
    return ProbabilistHermite(n).polynomial(x)
end
HermiteHe(x; n) = HermiteHe(n, x)

"""
    He(n, x)
    He(x; n)

Evaluates the nth probabilist's Hermite polynomial at x ∈ ℝ.
"""
He(n, x) = HermiteHe(n, x)
He(x; n) = HermiteHe(x; n)

# Functions for getting properties 
"""
    coefficients(He)

Determines the coefficients of a given probabilist's Hermite polynomial.
"""
coefficients(He::ProbabilistHermite) = He.coefficients 

"""
    weight(He)

Determines the weight function of a given probabilist's Hermite polynomial.
"""
weight(He::ProbabilistHermite) = He.weight 

"""
    interval(He)

Determines the interval of orthogonality of a given probabilsit's Hermite
polynomial.
"""
interval(He::ProbabilistHermite) = He.interval 

"""
    innerproduct(He1, He2)

Calculates the inner product of the probabilist's Hermite polynomials
`He1` and `He2` using the orthogonality relation 

``\\int_{-\\inf}^{\\inf} He_m(x) He_n(x) e^{-(x/2)^2} \\: \\mathrm{d}x = \\sqrt{2\\pi} n! \\delta_{nm}``.
"""
function innerproduct(He1::ProbabilistHermite, He2::ProbabilistHermite)
    n, m = He1.n, He2.m
    if n != m
        return 0
    else
        if n ≤ 20
            n! = factorial(n)
        else
            n! = gamma(n + 1)
        end
        return √(2π) * n!
    end
end

# ----------------------------------------------------------------------------------------------------------------------
# Exports 
# ----------------------------------------------------------------------------------------------------------------------

export Hermite 
export HermiteH 
export H 
export ProbabilistHermite 
export HermiteHe 
export He 
export coefficients 
export weight 
export interval 
export innerproduct

# ----------------------------------------------------------------------------------------------------------------------