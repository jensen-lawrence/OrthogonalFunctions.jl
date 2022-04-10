# ----------------------------------------------------------------------------------------------------------------------
# Hermite Polynomials 
# ----------------------------------------------------------------------------------------------------------------------

"""
    hermite_coefficients_direct(n)

Calculates the coefficients of the nth Hermite polynomial using the power series
representation of Hermite polynomials,

``H_n(x) = n! \\sum_{k=0}^{\\lfloor n/2 \\rfloor} \\frac{(-1)^k}{k!(n - 2k)!} (2x)^{n - 2k}``.
"""
function hermite_coefficients_direct(n::Int)
    nonzero = factorial(n) .* [((-1.0)^k * 2.0^(n - 2k))/(factorial(k) * factorial(n - 2k)) for k ∈ floor(n/2):-1:0]
    nonzero = map(x -> round(Int, x), nonzero)
    coefficients = zeros(Int, n + 1)
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
    hermite_coefficients_recursive(n)

Calculates the coefficients of the nth Hermite polynomial using the Hermite
polynomial recurrence relation

``H_n(x) = 2xH_{n-1}(x) - 2(n - 1)H_{n-2}(x)``.
"""
function hermite_coefficients_recursive(n::Int)
    last1 = zeros(n + 1)
    last2 = zeros(n + 1)
    last1[2:end] = hermite_coefficients(n - 1)
    last2[1:end-2] = hermite_coefficients(n - 2)
    coefficients = (2 * last1) - (2 * (n - 1) * last2)
    if maximum(coefficients) ≤ typemax(Int)
        coefficients = map(x -> round(Int, x), coefficients)
    end
    return SVector{n+1}(coefficients)
end

"""
    hermite_coefficients(n)

Calculates the coefficients of the nth Hermite polynomial. Switches from
`hermite_coefficients_direct` to `hermite_coefficients_recursive` for n > 20
due to integer overflow. Coefficients are ordered from lowest to highest
degree for compatibility with `evalpoly`.
"""
function hermite_coefficients(n::Int)
    if n ≤ 20
        return hermite_coefficients_direct(n)
    else
        return hermite_coefficients_recursive(n)
    end
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
    probabilist_hermite_coefficients_direct(n)

Calculates the coefficients of the nth probabilist's Hermite polynomial using
the power series representation of probabilist's Hermite polynomials,

``He_n(x) = n! \\sum_{k=0}^{\\lfloor n/2 \\rfloor} \\frac{(-1)^k}{k!(n - 2k)!} \\frac{x^{n - 2k}}{2^k}``.
"""
function probabilist_hermite_coefficients_direct(n::Int)
    nonzero = factorial(n) .* [((-1.0)^k)/(2.0^k * factorial(k) * factorial(n - 2k)) for k ∈ floor(n/2):-1:0]
    nonzero = map(x -> round(Int, x), nonzero)
    coefficients = zeros(Int, n + 1)
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
    probabilist_hermite_coefficients_recursive(n)

Calculates the coefficients of the nth probabilist's Hermite polynomial using
the probabilist's Hermite polynomial recurrence relation

``He_n(x) = xHe_{n-1}(x) - (n - 1)He_{n-2}(x)``.
"""
function probabilist_hermite_coefficients_recursive(n::Int)
    last1 = zeros(n + 1)
    last2 = zeros(n + 1)
    last1[2:end] = probabilist_hermite_coefficients(n - 1)
    last2[1:end-2] = probabilist_hermite_coefficients(n - 2)
    coefficients = last1 - ((n - 1) * last2)
    if maximum(coefficients) ≤ typemax(Int)
        coefficients = map(x -> round(Int, x), coefficients)
    end
    return SVector{n+1}(coefficients)
end

"""
    probabilist_hermite_coefficients(n)

Calculates the coefficients of the nth probabilist's Hermite polynomial.
Switches from `probabilist_hermite_coefficients_direct` to
`probabilist_hermite_coefficients_recursive` for n > 20 due to integer
overflow. Coefficients are ordered from lowest to highest degree for
compatibility with `evalpoly`.
"""
function probabilist_hermite_coefficients(n::Int)
    if n ≤ 20
        return probabilist_hermite_coefficients_direct(n)
    else
        return probabilist_hermite_coefficients_recursive(n)
    end
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