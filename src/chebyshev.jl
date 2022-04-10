# ----------------------------------------------------------------------------------------------------------------------
# Chebyshev Polynomials of the First Kind
# ----------------------------------------------------------------------------------------------------------------------

"""
    chebyshev_first_coefficients(n)

Calculates the coefficients of the nth Chebyshev polynomial of the first kind 
using the power series representation of Chebyshev polynomials of the first kind,

``T_n(x) = \\frac{n}{2} \\sum_{k=0}^{\\lfloor n/2 \\rfloor} \\frac{(-1)^k(n - k - 1)!}{k!(n - 2k)!} (2x)^{n - 2k}``.

Coefficients are ordered from lowest to highest degree for compatibility
with `evalpoly`.
"""
function chebyshev_first_coefficients(n::Int)
    if 1 ≤ n ≤ 20
        nonzero = (n/2) .* [((-1.0)^k * 2.0^(n - 2k) * factorial(n - k - 1))/(factorial(k) * factorial(n - 2k)) for k ∈ floor(n/2):-1:0]
    elseif n == 0
        nonzero = [0]
    else
        nonzero = (n/2) .* [((-1.0)^k * 2.0^(n - 2k) * gamma(n - k))/(gamma(k + 1) * gamma(n - 2k + 1)) for k ∈ floor(n/2):-1:0]
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
    ChebyshevFirst(n)

Representation of the nth Chebyshev polynomial of the first kind. Available
methods include:
- `n`: the polynomial degree 
- `coefficients`: the polynomial coefficients 
- `polynomial`: the polynomial function 
- `weight`: the Chebyshev polynomial of the first kind weight function, ``w(x) = \\frac{1}{\\sqrt{1 - x^2}}``
- `interval`: the Chebyshev polynomial of the first kind orthogonality interval, ``[-1, 1]``
"""
struct ChebyshevFirst <: AbstractChebyshev
    n::Int
    coefficients
    polynomial
    weight
    interval

    function ChebyshevFirst(n)
        coefficients = chebyshev_first_coefficients(n)
        polynomial = x -> evalpoly(x, coefficients)
        weight = x -> 1/√(1 - x^2)
        interval = (-1, 1)
        new(n, coefficients, polynomial, weight, interval)
    end
end

"""
    ChebyshevT(n, x)
    ChebyshevT(x; n)

Evaluates the nth Chebyshev polynomial of the first kind at x ∈ ℝ.
"""
function ChebyshevT(n::Int, x::T) where {T<:Real}
    return ChebyshevFirst(n).polynomial(x)
end
ChebyshevT(x; n) = ChebyshevT(n, x)

"""
    T(n, x)
    T(x; n)

Evaluates the nth Chebyshev polynomial of the first kind at x ∈ ℝ.
"""
T(n, x) = ChebyshevT(n, x)
T(x; n) = ChebyshevT(x; n)

"""
    coefficients(T)

Determines the coefficients of a given Chebyshev polynomial of the first kind.
"""
coefficients(T::ChebyshevFirst) = T.coefficients 

"""
    weight(T)

Determines the weight function of a given Chebyshev polynomial of the first kind.
"""
weight(T::ChebyshevFirst) = T.weight 

"""
    interval(T)

Determines the interval of orthogonality of a given Chebyshev polynomial of the first kind.
"""
interval(T::ChebyshevFirst) = T.interval

"""
    innerproduct(T1, T2)

Calculates the inner product of the Chebyshev polynomials of the first kind
`T1` and `T2` using the orthogonality relation 

``\\int_{-1}^{1} T_m(x) T_n(x) \\frac{1}{\\sqrt{1 - x^2}} \\: \\mathrm{d}x = \\frac{\\pi}{2}}\\delta_{nm}``,

or ``\\pi`` if ``n = m = 0``.
"""
function innerproduct(T1::ChebyshevFirst, T2::ChebyshevFirst)
    n, m = T1.n ,T2.n 
    if n != m 
        return 0
    elseif n == m == 0
        return π 
    else 
        return π/2
    end
end

# ----------------------------------------------------------------------------------------------------------------------
# Chebyshev Polynomials of the Second Kind
# ----------------------------------------------------------------------------------------------------------------------

"""
    chebyshev_second_coefficients(n)

Calculates the coefficients of the nth Chebyshev polynomial of the second kind 
using the power series representation of Chebyshev polynomials of the second kind,

``U_n(x) = \\sum_{k=0}^{\\lfloor n/2 \\rfloor} (-1)^k \\binom{n - k, k} (2x)^{n - 2k}``.

Coefficients are ordered from lowest to highest degree for compatibility
with `evalpoly`.
"""
function chebyshev_second_coefficients(n::Int)
    nonzero = [(-1.0)^k * 2.0^(n - 2k) * binomial(Int(n - k), Int(k)) for k ∈ floor(n/2):-1:0]
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
    ChebyshevSecond(n)

Representation of the nth Chebyshev polynomial of the second kind. Available
methods include:
- `n`: the polynomial degree 
- `coefficients`: the polynomial coefficients 
- `polynomial`: the polynomial function 
- `weight`: the Chebyshev polynomial of the second kind weight function, ``w(x) = \\sqrt{1 - x^2}``
- `interval`: the Chebyshev polynomial of the second kind orthogonality interval, ``[-1, 1]``
"""
struct ChebyshevSecond <: AbstractChebyshev
    n::Int
    coefficients
    polynomial
    weight
    interval

    function ChebyshevSecond(n)
        coefficients = chebyshev_second_coefficients(n)
        polynomial = x -> evalpoly(x, coefficients)
        weight = x -> (1 - x^2)
        interval = (-1, 1)
        new(n, coefficients, polynomial, weight, interval)
    end
end

"""
    ChebyshevU(n, x)
    ChebyshevU(x; n)

Evaluates the nth Chebyshev polynomial of the second kind at x ∈ ℝ.
"""
function ChebyshevU(n::Int, x::T) where {T<:Real}
    return ChebyshevSecond(n).polynomial(x)
end
ChebyshevU(x; n) = ChebyshevU(n, x)

"""
    T(n, x)
    T(x; n)

Evaluates the nth Chebyshev polynomial of the second kind at x ∈ ℝ.
"""
U(n, x) = ChebyshevU(n, x)
U(x; n) = ChebyshevU(x; n)

"""
    coefficients(U)

Determines the coefficients of a given Chebyshev polynomial of the second kind.
"""
coefficients(U::ChebyshevSecond) = U.coefficients 

"""
    weight(U)

Determines the weight function of a given Chebyshev polynomial of the second kind.
"""
weight(U::ChebyshevSecond) = U.weight 

"""
    interval(U)

Determines the interval of orthogonality of a given Chebyshev polynomial of the second kind.
"""
interval(U::ChebyshevSecond) = U.interval

"""
    innerproduct(U1, U2)

Calculates the inner product of the Chebyshev polynomials of the second kind
`U1` and `U2` using the orthogonality relation 

``\\int_{-1}^{1} U_m(x) U_n(x) \\sqrt{1 - x^2} \\: \\mathrm{d}x = \\frac{\\pi}{2}\\delta_{nm}``.
"""
function innerproduct(U1::ChebyshevSecond, U2::ChebyshevSecond)
    n, m = U1.n ,U2.n 
    if n != m 
        return 0
    else 
        return π/2
    end
end

# ----------------------------------------------------------------------------------------------------------------------
# Exports
# ----------------------------------------------------------------------------------------------------------------------

export ChebyshevFirst 
export ChebyshevT 
export T 
export ChebyshevSecond 
export ChebyshevU 
export U 
export coefficients 
export weight 
export interval 
export innerproduct

# ----------------------------------------------------------------------------------------------------------------------