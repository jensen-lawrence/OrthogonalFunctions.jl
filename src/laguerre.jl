# ----------------------------------------------------------------------------------------------------------------------
# Laguerre Polynomials
# ----------------------------------------------------------------------------------------------------------------------

"""
    laguerre_coefficients(n)

Calculates the coefficients of the nth Laguerre polynomial using the power
series representation of Laguerre polynomials,

``L_n(x) = \\sum_{k=0}^n (-1)^k \\binom{n, k} \frac{x^k}{k!}``

Coefficients are ordered from lowest to highest degree for compatibility
with `evalpoly`.
"""
function laguerre_coefficients(n::Int)
    coefficients = zeros(n + 1)
    for k ∈ 0:n
        if k ≤ 20
            coefficients[k+1] = ((-1)^k)/factorial(k) * binomial(n, k)
        else
            coefficients[k+1] = ((-1)^k)/gamma(k + 1) * binomial(n, k)
        end
    end
    return SVector{n+1}(coefficients)
end

"""
    Laguerre(n)

Representation of the nth Laguerre polynomial. Available methods include:
- `n`: the polynomial degree 
- `coefficients`: the polynomial coefficients 
- `polynomial`: the polynomial function 
- `weight`: the associated Laguerre polynomial weight function, ``w(x) = e^{-x}``
- `interval`: the associated Laguerre polynomial orthogonality interval, ``[0, \\infty)``
"""
struct Laguerre <: AbstractLaguerre
    n::Int 
    coefficients
    polynomial
    weight
    interval

    function Laguerre(n)
        coefficients = laguerre_coefficients(n)
        polynomial = x -> evalpoly(x, coefficients)
        weight = x -> exp(-x)
        interval = (0, Inf)
        new(n, coefficients, polynomial, weight, interval)
    end
end

"""
    LaguerreL(n, x)
    LaguerreL(x; n)

Evaluates the nth Laguerre polynomial at x ∈ ℝ.
"""
function LaguerreL(n::Int, x::T) where {T<:Real}
    return Laguerre(n).polynomial(x)
end
LaguerreL(x; n) = LaguerreL(n, x)

"""
    L(n, x)
    L(x; n)

 Evaluates the nth Laguerre polynomial at x ∈ ℝ.
"""
L(n, x) = LaguerreL(n, x)
L(x; n) = LaguerreL(x; n)

"""
    coefficients(L)

Determines the coefficients of a given Laguerre polynomial.
"""
coefficients(L::Laguerre) = L.coefficients 

"""
    weight(L)

Determines the weight function of a given Laguerre polynomial.
"""
weight(L::Laguerre) = L.weight 

"""
    interval(L)

Determines the interval of orthogonality of a given Laguerre polynomial.
"""
interval(L::Laguerre) = L.interval

"""
    innerproduct(L1, L2)

Calculates the inner product of the Laguerre polynomials `L1` and `L2`
using the orthogonality relation

``\\int_0^{\\inftyty} L_m^(x) L_n^(x) e^{-x} \\: \\mathrm{d}x = \\delta_{nm}``
"""
function innerproduct(L1::Laguerre, L2::Laguerre)
    n, m = L1.n, L2.n 
    if n != m 
        return 0
    else
        return 1
    end
end

# ----------------------------------------------------------------------------------------------------------------------
# Associated Laguerre Polynomials
# ----------------------------------------------------------------------------------------------------------------------

"""
    associated_laguerre_coefficients(n, α)

Calculates the coefficients of the (n, α) associated Laguerre polynomial using
the power series representation of associated Laguerre polynomials,

``L_n^{(\\alpha)}(x) = \\sum_{k=0}^n (-1)^k \\binom{n + α, n - k} \frac{x^k}{k!}``

Coefficients are ordered from lowest to highest degree for compatibility
with `evalpoly`.
"""
function associated_laguerre_coefficients(n::Int, α::T) where {T<:Real}
    coefficients = zeros(n + 1)
    for k ∈ 0:n
        if k ≤ 20
            coefficients[k+1] = ((-1)^k)/factorial(k) * binomial(n + α, n - k)
        else
            coefficients[k+1] = ((-1)^k)/gamma(k + 1) * binomial(n + α, n - k)
        end
    end
    return SVector{n+1}(coefficients)
end

"""
    AssociatedLaguerre(n, α)

Representation of the (n, α) associated Laguerre polynomial. Available methods
include:
- `n`: the polynomial degree
- `α`: the associated parameter 
- `coefficients`: the polynomial coefficients 
- `polynomial`: the polynomial function 
- `weight`: the associated Laguerre polynomial weight function, ``w(x) = x^{\\alpha}e^{-x}``
- `interval`: the associated Laguerre polynomial orthogonality interval, ``[0, \\infty)``

Note that `AssociatedLaguerre(n, 0)` = `Laguerre(n)`.
"""
struct AssociatedLaguerre{T<:Real} <: AbstractLaguerre
    n::Int 
    α::T 
    coefficients
    polynomial
    weight
    interval

    function AssociatedLaguerre(n, α)
        coefficients = associated_laguerre_coefficients(n, α)
        polynomial = x -> evalpoly(x, coefficients)
        weight = x -> (x^α)*exp(-x)
        interval = (0, Inf)
        new{T}(n, α, coefficients, polynomial, weight, interval)
    end
end

"""
    LaguerreLα(n, α, x)
    LaguerreLα(x; n, α)

Evaluates the (n, α) associated Laguerre polynomial at x ∈ ℝ.
"""
function LaguerreLα(n::Int, α::T1, x::T2) where {T1,T2<:Real}
    return AssociatedLaguerre(n, α).polynomial(x)
end
LaguerreLα(x; n, α) = LaguerreLα(n, α, x)

"""
    Lα(n, α, x)
    Lα(x; n, α)

 Evaluates the (n, α) associated Laguerre polynomial at x ∈ ℝ.
"""
Lα(n, α, x) = LaguerreLα(n, α, x)
Lα(x; n, α) = LaguerreLα(x; n, α)

"""
    coefficients(L)

Determines the coefficients of a given associated Laguerre polynomial.
"""
coefficients(L::AssociatedLaguerre) = L.coefficients 

"""
    weight(L)

Determines the weight function of a given associated Laguerre polynomial.
"""
weight(L::AssociatedLaguerre) = L.weight 

"""
    interval(L)

Determines the interval of orthogonality of a given associated Laguerre polynomial.
"""
interval(L::AssociatedLaguerre) = L.interval

"""
    innerproduct(L1, L2)

Calculates the inner product of the associated Laguerre polynomials `L1` and `L2`
using the orthogonality relation

``\\int_0^{\\inftyty} L_m^{(\\alpha)}(x) L_n^{(\\alpha)}(x) x^{\\alpha}e^{-x} \\: \\mathrm{d}x = \\frac{\\Gamma(n + α + 1)}{n!}\\delta_{nm}``
"""
function innerproduct(L1::AssociatedLaguerre, L2::AssociatedLaguerre)
    @assert L1.α == L2.α "Associated Laguerre polynomials must have the same α."
    n, m = L1.n, L2.n 
    if n != m 
        return 0
    else
        if n ≤ 20
            n! = factorial(n)
        else
            n! = gamma(n + 1)
        end
        return gamma(n + α + 1)/n! 
    end
end

# ----------------------------------------------------------------------------------------------------------------------
# Exports 
# ----------------------------------------------------------------------------------------------------------------------

export Laguerre 
export LaguerreL 
export L 
export AssociatedLaguerre 
export LaguerreLα
export Lα
export coefficients 
export weight 
export interval 
export innerproduct

# ----------------------------------------------------------------------------------------------------------------------