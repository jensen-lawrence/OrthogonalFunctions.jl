# ----------------------------------------------------------------------------------------------------------------------
# Hermite Polynomial
# ----------------------------------------------------------------------------------------------------------------------

"""
    HermitePolynomial(n)

Representation of the nth Hermite polynomial, where n ∈ ℤ⁺.

Available methods:
- `HermitePolynomial(n).n` returns n.
- `HermitePolynomial(n).f` returns the polyomial function.

`HermitePolynomial(n)` is directly callable as `HermitePolynomial(n)(x)`:
- if `typeof(x) <: Real`, this evaluates to the nth Hermite polynomial
evaluated at x.
- if `typeof(x) = SymPy.Sym`, this evaluates to the nth Hermite polynomial as
a symbolic expression with x as the variable.
- If `typeof(x) <: AbstractString`, this evaluates to the nth Hermite polynomial
as a LaTeX string with x as the variable.

Aliases: `HermiteH`, `H`.
"""
struct HermitePolynomial <: AbstractHermitePolynomial
    n::Integer
    f::Function

    function HermitePolynomial(n)
        n ≥ 0 || error("n must be non-negative.")
        f = x -> hermite_polynomial(n, x)
        new(n, f)
    end
end

# Make HermitePolynomial(n) callable 
function (H::HermitePolynomial)(x::T) where {T<:Union{Real,SymPy.Sym,AbstractString}}
    return H.f(x)
end

# HermitePolynomial(n) equality 
import Base.:(==)
==(H1::HermitePolynomial, H2::HermitePolynomial) = H1.n == H2.n

"""
    HermiteH(n)

Alias of `HermitePolynomial(n)`.
"""
function HermiteH(n::Integer)
    return HermitePolynomial(n)
end

"""
    H(n)

Alias of `HermitePolynomial(n)`.
"""
function H(n::Integer)
    return HermitePolynomial(n)
end

# ----------------------------------------------------------------------------------------------------------------------
# Differentiated Hermite Polynomial
# ----------------------------------------------------------------------------------------------------------------------

"""
    DifferentiatedHermitePolynomial(n, k)

Representation of the nth Hermite polynomial differentiated k times,
where n, k ∈ ℤ⁺.

Available methods:
- `DifferentiatedHermitePolynomial(n, k).n` returns n.
- `DifferentiatedHermitePolynomial(n, k).n` returns k.
- `DifferentiatedHermitePolynomial(n, k).f` returns the polyomial function.

`DifferentiatedHermitePolynomial(n, k)` is directly callable as
`DifferentiatedHermitePolynomial(n)(x, k)`:
- if `typeof(x) <: Real`, this evaluates to the kth derivative of the nth
Hermite polynomial evaluated at x.
- if `typeof(x) = SymPy.Sym`, this evaluates to the kth derivative of the nth
Hermite polynomial as a symbolic expression with x as the variable.
- If `typeof(x) <: AbstractString`, this evaluates to the kth derivative of the
nth Hermite polynomial as a LaTeX string with x as the variable.

Aliases: `dHermiteH`, `dH`.
"""
struct DifferentiatedHermitePolynomial
    n::Integer
    k::Integer
    f::Function

    function DifferentiatedHermitePolynomial(n, k)
        n ≥ 0 || error("n must be non-negative.")
        k ≥ 0 || error("k must be non-negative.")
        f = x -> differentiated_hermite_polynomial(n, k, x)
        new(n, k, f)
    end
end

# Make DifferentiatedHermitePolynomial(n) callable 
function (dH::DifferentiatedHermitePolynomial)(x::T) where {T<:Union{Real,SymPy.Sym,AbstractString}}
    return dH.f(x)
end

# DifferentiatedHermitePolynomial(n) equality 
==(H1::DifferentiatedHermitePolynomial, H2::DifferentiatedHermitePolynomial) = (H1.n == H2.n) && (H1.k == H2.k)

"""
    dHermiteH(n, k)

Alias of `DifferentiatedHermitePolynomial(n)`.
"""
function dHermiteH(n::Integer, k::Integer)
    return DifferentiatedHermitePolynomial(n, k)
end

"""
    dH(n, k)

Alias of `DifferentiatedHermitePolynomial(n)`.
"""
function dH(n::Integer, k::Integer)
    return DifferentiatedHermitePolynomial(n, k)
end

# ----------------------------------------------------------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------------------------------------------------------

"""
    hermite_polynomial(n, x)

If `typeof(x) <: Real`, returns the nth Hermite polynomial evaluated at x.

If `typeof(x) = SymPy.Sym`, returns the nth Hermite polynomial as a symbolic
expression with x as the variable.

If `typeof(x) <: AbstractString`, returns the nth Hermite polynomial as a 
LaTeX string with x as the variable.
"""
function hermite_polynomial(n::Integer, x::T) where {T<:Real}
    Σ = zero(x)
    @simd for k ∈ 0:fld(n, 2)
        Σ += ((-1)^k)/(fac(k) * fac(n - 2k)) * (2x)^(n - 2k)
    end
    return fac(n) * Σ
end

function hermite_polynomial(n::Integer, x::SymPy.Sym)
    coeff = Int((-1.0)^n) * exp(x^2)
    expr = coeff * diff(exp(-x^2), x, n)
    return factor(expr)
end

function hermite_polynomial(n::Integer, x::AbstractString)
    sym = symbols(x)
    expr = latexify(hermite_polynomial(n, sym))
    expr = replace(expr, "\\cdot " => "")
    return latexstring(expr[1], "H_{$n}($x) = ", expr[2:end])
end

"""
    differentiated_hermite_polynomial(n, x)

If `typeof(x) <: Real`, returns the nth Hermite polynomial differentiated
k times evaluated at x.

If `typeof(x) = SymPy.Sym`, returns the nth Hermite polynomial differentiated
k times as a symbolic expression with x as the variable.

If `typeof(x) <: AbstractString`, returns the nth Hermite polynomial
differentiated k times as a LaTeX string with x as the variable.
"""
function differentiated_hermite_polynomial(n::Integer, k::Integer, x::T) where {T<:Real}
    if n ≥ k 
        Hₙ₋ₖ = hermite_polynomial(n - k, x)
        return 2^k * fac(n)/fac(n - k) * Hₙ₋ₖ
    else
        return zero(x)
    end
end

function differentiated_hermite_polynomial(n::Integer, k::Integer, x::SymPy.Sym)
    coeff = Int((-1.0)^(n - k) * 2^k * fac(n)/fac(n - k)) * exp(x^2)
    expr = coeff * diff(exp(-x^2), x, n - k)
    return factor(expr)
end

function differentiated_hermite_polynomial(n::Integer, k::Integer, x::AbstractString)
    sym = symbols(x)
    expr = latexify(differentiated_hermite_polynomial(n, k, sym))
    expr = replace(expr, "\\cdot " => "")
    return latexstring(expr[1], "\\frac{\\mathrm{d}^{$k}H_{$n}}{\\mathrm{d}x^{$k}} = ", expr[2:end])
end

"""
    hermite_polynomial_weight(x)

If `typeof(x) <: Real`, returns the orthogonality weight function for Hermite 
polynomials evaluated at x.

If `typeof(x) = SymPy.Sym`, returns the orthogonality weight function for
Hermite polynomials as a symbolic expression with x as the variable.

If `typeof(x) <: AbstractString`, returns the orthogonality weight function for
Hermite polynomials as a LaTeX string with x as the variable.
"""
function hermite_polynomial_weight(x::T) where {T<:Real}
    return exp(-x^2)
end

function hermite_polynomial_weight(x::SymPy.Sym)
    return exp(-x^2)
end

function hermite_polynomial_weight(x::AbstractString)
    return L"w(%$x) = \exp{-%$x^2}"
end

"""
    weight(H)

Returns the orthogonality weight function for Hermite polynomials.
"""
function weight(H::HermitePolynomial)
    return hermite_polynomial_weight
end

"""
    innerproduct(H1, H2)

Returns the inner product of the Hermite polynomials H1 and H2 using the
orthogonality relation for Hermite polynomials.
"""
function innerproduct(H1::HermitePolynomial, H2::HermitePolynomial)
    n, m = H1.n, H2.n 
    return √π * 2^n * fac(n) * δ(n, m)
end

"""
    derivative(H, k)

Returns the kth derivative of the Hermite polynomial H.
"""
function derivative(H::HermitePolynomial, k::Integer)
    return DifferentiatedHermitePolynomial(H.n, k)
end

# ----------------------------------------------------------------------------------------------------------------------