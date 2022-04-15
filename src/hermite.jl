# ----------------------------------------------------------------------------------------------------------------------
# Hermite Polynomials
# ----------------------------------------------------------------------------------------------------------------------

"""
    hermite_H_f(n, x)

Evaluates the nth Hermite polynomial at x ∈ ℝ.
"""
function hermite_H_f(n::Int, x::R) where {R<:Real}
    H = sum([((-1)^k)/(genfac(k) * genfac(n - 2k)) * (2x)^(n - 2k) for k ∈ 0:floor(Int, n/2)])
    H *= genfac(n)
    return H 
end

"""
    hermite_H_w(x)

Evaluates the orthogonality weight function of the nth Hermite polynomial
at x ∈ ℝ.
"""
function hermite_H_w(x::R) where {R<:Real}
    return exp(-x^2)
end

"""
    HermiteH(n)

Representation of the nth Hermite polynomial, Hₙ, where n ∈ ℤ⁺.

Available methods:
- `HermiteH(n).n` returns n (the polynomial degree).
- `HermiteH(n).f` returns the polyomial function.
- `HermiteH(n).w` returns the orthogonality weight function.
- `HermiteH(n).I` returns the orthogonality interval.

`HermiteH(n)` is directly callable, so Hₙ(x) can be evaluated using 
`HermiteH(n)(x)`. Hₙ(x) can also be evaluated using `HermiteH(n).f(x)`.
Similarly, the orthogonality weight function w(x) can be evaluated using 
`HermiteH(n).w(x)`.

Aliases: `H`, `Hermite`.
"""
struct HermiteH <: AbstractHermite 
    n::Int 
    f 
    w 
    I 

    function HermiteH(n)
        n ≥ 0 || error("n must be non-negative.")
        f = x -> hermite_H_f(n, x)
        w = x -> hermite_H_w(x)
        I = (-Inf, Inf)
        new(n, f, w, I)
    end
end

# Make HermiteH(n) callable 
function (H::HermiteH)(x::R) where {R<:Real}
    return H.f(x)
end

"""
    H(n)

Alias of `HermiteH`.
"""
function H(n::Int)
    return HermiteH(n)
end

"""
    Hermite(n)

Alias of `HermiteH`.
"""
function Hermite(n::Int)
    return HermiteH(n)
end

"""
    innerproduct(H1, H2)

Returns the inner product of `H1` and `H2` using the orthogonality relation 
for Hermite polynomials.
"""
function innerproduct(H1::HermiteH, H2::HermiteH)
    n, m = H1.n, H2.n 
    return √π * 2^n * genfac(n) * δ(n, m)
end

# ----------------------------------------------------------------------------------------------------------------------
# Probabilist's Hermite Polynomials
# ----------------------------------------------------------------------------------------------------------------------

"""
    hermite_He_f(n, x)

Evaluates the nth probabilist's Hermite polynomial at x ∈ ℝ.
"""
function hermite_He_f(n::Int, x::R) where {R<:Real}
    He = sum([((-1)^k)/(genfac(k) * genfac(n - 2k)) * ((2x)^(n - 2k))/(2^k) for k ∈ 0:floor(Int, n/2)])
    He *= genfac(n)
    return He 
end

"""
    hermite_He_w(x)

Evaluates the orthogonality weight function of the nth probabilist's
Hermite polynomial at x ∈ ℝ.
"""
function hermite_He_w(x::R) where {R<:Real}
    return exp(-(x^2)/2)
end

"""
    HermiteHe(n)

Representation of the nth Hermite polynomial, Heₙ, where n ∈ ℤ⁺.

Available methods:
- `HermiteHe(n).n` returns n (the polynomial degree).
- `HermiteHe(n).f` returns the polyomial function.
- `HermiteHe(n).w` returns the orthogonality weight function.
- `HermiteHe(n).I` returns the orthogonality interval.

`HermiteHe(n)` is directly callable, so Heₙ(x) can be evaluated using 
`HermiteHe(n)(x)`. Heₙ(x) can also be evaluated using `HermiteHe(n).f(x)`.
Similarly, the orthogonality weight function w(x) can be evaluated using 
`HermiteHe(n).w(x)`.

Aliases: `He`, `ProbabilistHermite`.
"""
struct HermiteHe <: AbstractHermite 
    n::Int 
    f 
    w 
    I 

    function HermiteHe(n)
        n ≥ 0 || error("n must be non-negative.")
        f = x -> hermite_He_f(n, x)
        w = x -> hermite_He_w(x)
        I = (-Inf, Inf)
        new(n, f, w, I)
    end
end

# Make HermiteHe(n) callable 
function (He::HermiteHe)(x::R) where {R<:Real}
    return He.f(x)
end

"""
    He(n)

Alias of `HermiteHe`.
"""
function He(n::Int)
    return HermiteHe(n)
end

"""
    ProbabilistHermite(n)

Alias of `HermiteHe`.
"""
function ProbabilistHermite(n::Int)
    return HermiteHe(n)
end

"""
    innerproduct(He1, He2)

Returns the inner product of `He1` and `He2` using the orthogonality relation 
for probabilist's Hermite polynomials.
"""
function innerproduct(He1::HermiteHe, He2::HermiteHe)
    n, m = He1.n, He2.n 
    return √(2π) * genfac(n) * δ(n, m)
end

# ----------------------------------------------------------------------------------------------------------------------
# General
# ----------------------------------------------------------------------------------------------------------------------

"""
    roots(F)

Returns the roots (zeros) of `F`.
"""
function roots(F::AbstractHermite)
    bound = -2*√(F.n)
    return find_zeros(F, (-bound, bound))
end

# ----------------------------------------------------------------------------------------------------------------------