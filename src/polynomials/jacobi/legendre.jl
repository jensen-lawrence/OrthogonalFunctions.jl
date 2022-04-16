# ----------------------------------------------------------------------------------------------------------------------
# Legendre Polynomials
# ----------------------------------------------------------------------------------------------------------------------

"""
    legendre_P_f(n, x)

Evaluates the nth Legendre polynomial at x ∈ ℝ.
"""
function legendre_P_f(n::Int, x::R) where {R<:Real}
    P = sum([genbinom(n, k)^2 * (x - 1)^(n - k) * (x + 1)^k for k ∈ 0:n])
    P *= 1/2^n 
    return P 
end

"""
    legendre_P_w(x)

Evaluates the orthogonality weight function of the nth Legendre polynomial 
at x ∈ ℝ.
"""
function legendre_P_w(x::R) where {R<:Real}
    return 1
end

"""
    LegendreP(n)

Representation of the nth Legendre polynomial, Tₙ, where n ∈ ℤ⁺.

Available methods:
- `LegendreP(n).n` returns n (the polynomial degree).
- `LegendreP(n).f` returns the polyomial function.
- `LegendreP(n).w` returns the orthogonality weight function.
- `LegendreP(n).I` returns the orthogonality interval.

`LegendreP(n)` is directly callable, so Pₙ(x) can be evaluated using 
`LegendreP(n)(x)`. Pₙ(x) can also be evaluated using `LegendreP(n).f(x)`.
Similarly, the orthogonality weight function w(x) can be evaluated using 
`LegendreP(n).w(x)`.

Aliases: `P`, `Legendre`.
"""
struct LegendreP <: AbstractJacobi
    n::Int 
    f
    w 
    I 

    function LegendreP(n)
        n ≥ 0 || error("n must be non-negative.")
        f = x -> legendre_P_f(n, x)
        w = x -> legendre_P_w(x)
        I = (-1, 1)
        new(n, f, w, I)
    end
end

# Make LegendreP(n) callable 
function (P::LegendreP)(x::R) where {R<:Real}
    return P.f(x)
end

"""
    P(n)

Alias of `LegendreP(n)`.
"""
function P(n::Int)
    return LegendreP(n)
end

"""
    Legendre(n)

Alias of `LegendreP(n)`.
"""
function Legendre(n::Int)
    return LegendreP(n)
end

"""
    innerproduct(P1, P2)

Returns the inner product of `P1` and `P2` using the orthogonality relation 
for Legendre polynomials.
"""
function innerproduct(P1::LegendreP, P2::LegendreP)
    n, m = P1.n, P2.n 
    return 2/(2n + 1) * δ(n, m)
end

# ----------------------------------------------------------------------------------------------------------------------