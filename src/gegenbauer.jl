# ----------------------------------------------------------------------------------------------------------------------
# Gegenbauer Polynomials
# ----------------------------------------------------------------------------------------------------------------------

"""
    gegenbauer_C_f(n, a, x)

Evaluates the (n, a) Gegenbauer polymomial at x ∈ ℝ.
"""
function gegenbauer_C_f(n::Int, a::R1, x::R2) where {R1,R2<:Real}
    C = sum([(-1)^k * (gamma(n - k + a)/(gamma(a) * genfac(k) * genfac(n - 2k))) * (2x)^(n - 2k) for k ∈ 0:floor(Int, n/2)])
    return C
end

"""
    gegenbauer_C_w(a, x)

Evaluates the orthogonality weight function of the (n, a) Gegenbauer polynomial
at x ∈ ℝ.
"""
function gegenbauer_C_w(a::R1, x::R2) where {R1,R2<:Real}
    return (1 - x^2)^(a - 1/2)
end

"""
    GegenbauerC(n, a)

Representation of the (n, a) Gegenbauer polynomial, Cₙᵃ, where n ∈ ℤ⁺ and
a ∈ ℝ > -1/2.

Available methods:
- `GegenbauerC(n, a).n` returns n (the polynomial degree).
- `GegenbauerC(n, a).a` returns a.
- `GegenbauerC(n, a).f` returns the polyomial function.
- `GegenbauerC(n, a).w` returns the orthogonality weight function.
- `GegenbauerC(n, a).I` returns the orthogonality interval.

`GegenbauerC(n, a)` is directly callable, so Cₙᵃ(x) can be evaluated using 
`GegenbauerC(n, a)(x)`. Cₙᵃ(x) can also be evaluated using`GegenbauerC(n, a).f(x)`.
Similarly, the orthogonality weight function w(x) can be evaluated using 
`GegenbauerC(n, a).w(x)`.

Aliases: `C`, `Gegenbauer`, `Ultraspherical`.
"""
struct GegenbauerC{R<:Real} <: OrthogonalPolynomial
    n::Int 
    a::R 
    f
    w 
    I 

    function GegenbauerC(n, a)
        n ≥ 0 || error("n must be non-negative.")
        a > -0.5 || error("a must be greater than -1/2.")
        f = x -> gegenbauer_C_f(n, a, x)
        w = x -> gegenbauer_C_w(a, x)
        I = (-1, 1)
        R = promote_type(typeof(a))
        new{R}(n, a, f, w, I)
    end
end

# Make GegenbauerC(n, a) callable 
function (C::GegenbauerC)(x::R) where {R<:Real}
    return C.f(x)
end

"""
    C(n, a)

Alias of `GegenbauerC(n, a)`.
"""
function C(n::Int, a::R) where {R<:Real}
    return GegenbauerC(n, a)
end

"""
    Gegenbauer(n, a)

Alias of `GegenbauerC(n, a)`.
"""
function Gegenbauer(n::Int, a::R) where {R<:Real}
    return GegenbauerC(n, a)
end

"""
    Ultraspherical(n, a)

Alias of `GegenbauerC(n, a)`.
"""
function Ultraspherical(n::Int, a::R) where {R<:Real}
    return GegenbauerC(n, a)
end

"""
    roots(C)

Returns the roots (zeros) of `C`.
"""
function roots(C::Gegenbauer)
    return find_zeros(C, (-1, 1))
end

"""
    innerproduct(C1, C2)

Returns the inner product of `C1` and `C2` using the orthogonality relation 
for Gegenbauer polynomials.
"""
function innerproduct(C1::GegenbauerC, C2::GegenbauerC)
    a1, a2 = C1.a, C2.a 
    a1 == a2 || error("C1 and C2 must have the same values of a.")
    a = a1 
    n, m = C1.n, C2.n 
    innerprod = (π * 2^(1 - 2a) * gamma(n + 2a))/(genfac(n) * (n + a) * gamma(a)^2) * δ(n, m)
    return innerprod
end

# ----------------------------------------------------------------------------------------------------------------------