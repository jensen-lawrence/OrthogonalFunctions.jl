# ----------------------------------------------------------------------------------------------------------------------
# Jacobi Polynomials
# ----------------------------------------------------------------------------------------------------------------------

"""
    jacobi_P_f(n, a, b, x)

Evaluates the (n, a, b) Jacobi polynomial at x ∈ ℝ.
"""
function jacobi_P_f(n::Int, a::R1, b::R2, x::R3) where {R1,R2,R3<:Real}
    P = sum([genbinom(n, k) * (gamma(a + b + n + k + 1)/gamma(a + k + 1)) * ((x - 1)/2)^k for k ∈ 0:n])
    P *= gamma(a + n + 1)/(genfac(n) * gamma(a + b + n + 1))
    return P
end

"""
    jacobi_P_w(a, b, x)

Evaluates the orthogonality weight function of the (n, a, b) Jacobi polynomial
at x ∈ ℝ.
"""
function jacobi_P_w(a::R1, b::R2, x::R3) where {R1,R2,R3<:Real}
    return (1 - x)^a * (1 + x)^b
end

"""
    JacobiP(n, a, b)

Representation of the (n, a, b) Jacobi polynomial, Pₙᵃᵇ, where n ∈ ℤ⁺ and 
a, b ∈ ℝ > -1.

Available methods:
- `JacobiP(n, a, b).n` returns n (the polynomial degree).
- `JacobiP(n, a, b).a` returns a.
- `JacobiP(n, a, b).b` returns b. 
- `JacobiP(n, a, b).f` returns the polyomial function.
- `JacobiP(n, a, b).w` returns the orthogonality weight function.
- `JacobiP(n, a, b).I` returns the orthogonality interval.

`JacobiP(n, a, b)` is directly callable, so Pₙᵃᵇ(x) can be evaluated using 
`JacobiP(n, a, b)(x)`. Pₙᵃᵇ(x) can also be evaluated using `JacobiP(n, a, b).f(x)`.
Similarly, the orthogonality weight function w(x) can be evaluated using 
`JacobiP(n, a, b).w(x)`.

Aliases: `P`, `Jacobi`.
"""
struct JacobiP{R1<:Real,R2<:Real} <: OrthogonalPolynomial 
    n::Int
    a::R1 
    b::R2 
    f
    w
    I 

    function JacobiP(n, a, b)
        n ≥ 0 || error("n must be non-negative.")
        a > -1 || error("a must be greater than -1.")
        b > -1 || error("b must be greater than -1.")
        f = x -> jacobi_P_f(n, a, b, x)
        w = x -> jacobi_P_w(a, b, x)
        I = (-1, 1)
        R1 = promote_type(typeof(a))
        R2 = promote_type(typeof(b))
        new{R1,R2}(n, a, b, f, w, I)
    end
end

# Make JacobiP(n, a, b) callable 
function (P::JacobiP)(x::R) where {R<:Real}
    return P.f(x)
end

"""
    P(n, a, b)

Alias of `JacobiP`.
"""
function P(n::Int, a::R1, b::R2) where {R1,R2<:Real}
    return JacobiP(n, a, b)
end

"""
    Jacobi(n, a, b)

Alias of `JacobiP`.
"""
function Jacobi(n::Int, a::R1, b::R2) where {R1,R2<:Real}
    return JacobiP(n, a, b)
end

"""
    degree(P)

Returns the degree of `P`.
"""
degree(P::JacobiP) = P.n 

"""
    polynomial(P)

Returns the polynomial function of `P`.
"""
polymomial(P::JacobiP) = P.f 

"""
    weight(P)

Returns the orthogonality weight function of `P`.
"""
weight(P::JacobiP) = P.w 

"""
    interval(P)

Returns the orthogonality interval of `P`.
"""
interval(P::JacobiP) = P.I 

"""
    roots(P)

Returns the roots (zeros) of `P`.
"""
function roots(P::JacobiP)
    return find_zeros(P, (-1, 1))
end

"""
    innerproduct(P1, P2)

Returns the inner product of `P1` and `P2` using the orthgonality relation 
for Jacobi polynomials.
"""
function innerproduct(P1::JacobiP, P2::JacobiP)
    a1, b1 = P1.a, P1.b 
    a2, b2 = P2.a, P2.b
    a1 == a2 || error("P1 and P2 must have the same values of a.")
    b1 == b2 || error("P1 and P2 must have the same values of b.")
    a, b = a1, b1 
    n, m = P1.n, P2.n 
    ξ = (2^(a + b + 1))/(2n + a + b + 1)
    χ = (gamma(n + a + 1) * gamma(n + b + 1))/(gamma(n + a + b + 1) * genfac(n))
    innerprod = ξ * χ * δ(n, m)
    return innerprod
end

"""
    raise(P)

Raises `P` from Pₙ to Pₙ₊₁.
"""
function raise(P::JacobiP)
    n, a, b = P.n, P.a, P.b 
    return JacobiP(n + 1, a, b)
end

"""
    lower(P)

Lowers `P` from Pₙ to Pₙ₋₁.
"""
function lower(P::JacobiP)
    n, a, b = P.n, P.a, P.b 
    n ≥ 1 || error("P₀ cannot be lowered further.")
    return JacobiP(n - 1, a, b)
end

# ----------------------------------------------------------------------------------------------------------------------