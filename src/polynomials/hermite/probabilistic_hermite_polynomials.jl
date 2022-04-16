# ----------------------------------------------------------------------------------------------------------------------
# Hermite Polynomial
# ----------------------------------------------------------------------------------------------------------------------

"""
    ProbabilisticHermitePolynomial(n)

Representation of the nth probabilistic Hermite polynomial, where n ∈ ℤ⁺.

Available methods:
- `ProbabilisticHermitePolynomial(n).n` returns n.
- `ProbabilisticHermitePolynomial(n).f` returns the polyomial function.

`ProbabilisticHermitePolynomial(n)` is directly callable as
`ProbabilisticHermitePolynomial(n)(x)`:
- if `typeof(x) <: Real`, this evaluates to the nth probabilistic Hermite
polynomial evaluated at x.
- if `typeof(x) = SymPy.Sym`, this evaluates to the nth probabilistic Hermite
polynomial as a symbolic expression with x as the variable.
- If `typeof(x) <: AbstractString`, this evaluates to the nth probabilistic
Hermite polynomial as a LaTeX string with x as the variable.

Aliases: `HermiteHe`, `He`.
"""
struct ProbabilisticHermitePolynomial <: AbstractHermitePolynomial
    n::Integer
    f::Function

    function ProbabilisticHermitePolynomial(n)
        n ≥ 0 || error("n must be non-negative.")
        f = x -> probabilistic_hermite_polynomial(n, x)
        new(n, f)
    end
end

# Make HermitePolynomial(n) callable 
function (He::ProbabilisticHermitePolynomial)(x::T) where {T<:Union{Real,SymPy.Sym,AbstractString}}
    return He.f(x)
end

# ProbabilisticHermitePolynomial(n) equality 
import Base.:(==)
==(He1::ProbabilisticHermitePolynomial, He2::ProbabilisticHermitePolynomial) = He1.n == He2.n

"""
    HermiteHe(n)

Alias of `ProbabilisticHermitePolynomial(n)`.
"""
function HermiteHe(n::Integer)
    return ProbabilisticHermitePolynomial(n)
end

"""
    He(n)

Alias of `ProbabilisticHermitePolynomial(n)`.
"""
function He(n::Integer)
    return ProbabilisticHermitePolynomial(n)
end

# ----------------------------------------------------------------------------------------------------------------------
# Differentiated Probabilistic Hermite Polynomial
# ----------------------------------------------------------------------------------------------------------------------

"""
    DifferentiatedProbabilisticHermitePolynomial(n, k)

Representation of the nth Hermite polynomial differentiated k times,
where n, k ∈ ℤ⁺.

Available methods:
- `DifferentiatedProbabilisticHermitePolynomial(n, k).n` returns n.
- `DifferentiatedProbabilisticHermitePolynomial(n, k).n` returns k.
- `DifferentiatedProbabilisticHermitePolynomial(n, k).f` returns the polyomial
function.

`DifferentiatedProbabilisticHermitePolynomial(n, k)` is directly callable as
`DifferentiatedProbabilisticHermitePolynomial(n)(x, k)`:
- if `typeof(x) <: Real`, this evaluates to the kth derivative of the nth
probabilistic Hermite polynomial evaluated at x.
- if `typeof(x) = SymPy.Sym`, this evaluates to the kth derivative of the nth
probabilistic Hermite polynomial as a symbolic expression with x as the variable.
- If `typeof(x) <: AbstractString`, this evaluates to the kth derivative of the
nth probabilistic Hermite polynomial as a LaTeX string with x as the variable.

Aliases: `dHermiteH`, `dH`.
"""
struct DifferentiatedProbabilisticHermitePolynomial
    n::Integer
    k::Integer
    f::Function

    function DifferentiatedProbabilisticHermitePolynomial(n, k)
        n ≥ 0 || error("n must be non-negative.")
        k ≥ 0 || error("k must be non-negative.")
        f = x -> differentiated_probabilistic_hermite_polynomial(n, k, x)
        new(n, k, f)
    end
end

# Make DifferentiatedProbabilisticHermitePolynomial(n) callable 
function (dHe::DifferentiatedProbabilisticHermitePolynomial)(x::T) where {T<:Union{Real,SymPy.Sym,AbstractString}}
    return dHe.f(x)
end

# DifferentiatedProbabilisticHermitePolynomial(n) equality 
==(He1::DifferentiatedProbabilisticHermitePolynomial, He2::DifferentiatedProbabilisticHermitePolynomial) = (He1.n == He2.n) && (He1.k == He2.k)

"""
    dHermiteHe(n, k)

Alias of `DifferentiatedProbabilisticHermitePolynomial(n)`.
"""
function dHermiteHe(n::Integer, k::Integer)
    return DifferentiatedProbabilisticHermitePolynomial(n, k)
end

"""
    dHe(n, k)

Alias of `DifferentiatedProbabilisticHermitePolynomial(n)`.
"""
function dHe(n::Integer, k::Integer)
    return DifferentiatedProbabilisticHermitePolynomial(n, k)
end

# ----------------------------------------------------------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------------------------------------------------------

"""
    probabilistic_hermite_polynomial(n, x)

If `typeof(x) <: Real`, returns the nth probabilistic Hermite polynomial
evaluated at x.

If `typeof(x) = SymPy.Sym`, returns the nth probabilistic Hermite polynomial
as a symbolic expression with x as the variable.

If `typeof(x) <: AbstractString`, returns the nth probabilistic Hermite
polynomial as a  LaTeX string with x as the variable.
"""
function probabilistic_hermite_polynomial(n::Integer, x::T) where {T<:Real}
    Σ = zero(x)
    @simd for k ∈ 0:fld(n, 2)
        Σ += ((-1)^k)/(fac(k) * fac(n - 2k) * 2^k) * x^(n - 2k)
    end
    return fac(n) * Σ
end

function probabilistic_hermite_polynomial(n::Integer, x::SymPy.Sym)
    coeff = Int((-1.0)^n) * exp((x^2)/2)
    expr = coeff * diff(exp(-(x^2)/2), x, n)
    return factor(expr)
end

function probabilistic_hermite_polynomial(n::Integer, x::AbstractString)
    sym = symbols(x)
    expr = latexify(probabilistic_hermite_polynomial(n, sym))
    expr = replace(expr, "\\cdot " => "")
    return latexstring(expr[1], "He_{$n}($x) = ", expr[2:end])
end

"""
    differentiated_probabilistic_hermite_polynomial(n, x)

If `typeof(x) <: Real`, returns the nth probabilistic Hermite polynomial
differentiated k times evaluated at x.

If `typeof(x) = SymPy.Sym`, returns the nth probabilistic Hermite polynomial
differentiated k times as a symbolic expression with x as the variable.

If `typeof(x) <: AbstractString`, returns the nth probabilistic Hermite
polynomial differentiated k times as a LaTeX string with x as the variable.
"""
function differentiated_probabilistic_hermite_polynomial(n::Integer, k::Integer, x::T) where {T<:Real}
    if n ≥ k 
        Hₙ₋ₖ = probabilistic_hermite_polynomial(n - k, x)
        return fac(n)/fac(n - k) * Hₙ₋ₖ
    else
        return zero(x)
    end
end

function differentiated_probabilistic_hermite_polynomial(n::Integer, k::Integer, x::SymPy.Sym)
    coeff = Int((-1.0)^(n - k) * fac(n)/fac(n - k)) * exp((x^2)/2)
    expr = coeff * diff(exp(-(x^2)/2), x, n - k)
    return factor(expr)
end

function differentiated_probabilistic_hermite_polynomial(n::Integer, k::Integer, x::AbstractString)
    sym = symbols(x)
    expr = latexify(differentiated_probabilistic_hermite_polynomial(n, k, sym))
    expr = replace(expr, "\\cdot " => "")
    return latexstring(expr[1], "\\frac{\\mathrm{d}^{$k}He_{$n}}{\\mathrm{d}x^{$k}} = ", expr[2:end])
end

"""
    probabilistic_hermite_polynomial_weight(x)

If `typeof(x) <: Real`, returns the orthogonality weight function for
probabilistic Hermite polynomials evaluated at x.

If `typeof(x) = SymPy.Sym`, returns the orthogonality weight function for
probabilistic Hermite polynomials as a symbolic expression with x as the
variable.

If `typeof(x) <: AbstractString`, returns the orthogonality weight function for
probabilistic Hermite polynomials as a LaTeX string with x as the variable.
"""
function probabilistic_hermite_polynomial_weight(x::T) where {T<:Real}
    return exp(-0.5 * (x^2))
end

function probabilistic_hermite_polynomial_weight(x::SymPy.Sym)
    return exp(-(x^2)/2)
end

function probabilistic_hermite_polynomial_weight(x::AbstractString)
    return L"w(%$x) = \exp{-(%$x^2)/2}"
end

"""
    weight(He)

Returns the orthogonality weight function for probabilistic Hermite
    polynomials.
"""
function weight(He::ProbabilisticHermitePolynomial)
    return probabilistic_hermite_polynomial_weight
end

"""
    innerproduct(He1, He2)

Returns the inner product of the Hermite polynomials He1 and He2 using the
orthogonality relation for probabilistic Hermite polynomials.
"""
function innerproduct(He1::ProbabilisticHermitePolynomial, He2::ProbabilisticHermitePolynomial)
    n, m = He1.n, He2.n 
    return √(2π) * fac(n) * δ(n, m)
end

"""
    derivative(He, k)

Returns the kth derivative of the probabilistic Hermite polynomial He.
"""
function derivative(He::ProbabilisticHermitePolynomial, k::Integer)
    return DifferentiatedProbabilisticHermitePolynomial(He.n, k)
end

# ----------------------------------------------------------------------------------------------------------------------