# ----------------------------------------------------------------------------------------------------------------------
# Laguerre Polynomials
# ----------------------------------------------------------------------------------------------------------------------

"""
    laguerre_L_f(n, x)

Evaluates the nth Laguerre polynomial at x ∈ ℝ.
"""
function laguerre_L_f(n::Int, x::R) where {R<:Real}
    L = sum([genbinom(n, k) * ((-1)^k)/genfac(k) * x^k for k ∈ 0:n])
    return L
end

"""
    laguerre_L_w(x)

Evaluates the orthogonality weight function of the nth Laguerre polynomial 
at x ∈ ℝ.
"""
function laguerre_L_w(x::R) where {R<:Real}
    return exp(-x)
end

"""
    LaguerreL(n)

Representation of the nth Laguerre polynomial, Lₙ, where n ∈ ℤ⁺.

Available methods:
- `LaguerreL(n).n` returns n (the polynomial degree).
- `LaguerreL(n).f` returns the polyomial function.
- `LaguerreL(n).w` returns the orthogonality weight function.
- `LaguerreL(n).I` returns the orthogonality interval.

`LaguerreL(n)` is directly callable, so Lₙ(x) can be evaluated using 
`LaguerreL(n)(x)`. Lₙ(x) can also be evaluated using `LaguerreL(n).f(x)`.
Similarly, the orthogonality weight function w(x) can be evaluated using 
`LaguerreL(n).w(x)`.

Aliases: `L`, `Laguerre`.
"""
struct LaguerreL <: AbstractLaguerre 
    n::Int 
    f 
    w 
    I 

    function LaguerreL(n)
        n ≥ 0 || error("n must be non-negative.")
        f = x -> laguerre_L_f(n, x)
        w = x -> laguerre_L_w(x)
        I = (0, Inf)
        new(n, f, w, I)
    end
end

# Make LaguerreL(n) callable 
function (L::LaguerreL)(x::R) where {R<:Real}
    return L.f(x)
end

"""
    L(n)

Alias of `LaguerreL`.
"""
function L(n::Int)
    return LaguerreL(n)
end

"""
    Laguerre(n)

Alias of `LaguerreL`.
"""
function Laguerre(n::Int)
    return LaguerreL(n)
end

"""
    innerproduct(L1, L2)

Returns the inner product of `L1` and `L2` using the orthogonality relation 
for Laguerre polynomials.
"""
function innerproduct(L1::LaguerreL, L2::LaguerreL)
    n, m = L1.n, L2.n 
    return δ(n, m)
end

# ----------------------------------------------------------------------------------------------------------------------
# General
# ----------------------------------------------------------------------------------------------------------------------

"""
    roots(F)

Returns the roots (zeros) of `F`.
"""
function roots(F::AbstractLaguerre)
    bound = 2*(2*F.n + 1)
    return find_zeros(F, (-bound, bound))
end

# ----------------------------------------------------------------------------------------------------------------------