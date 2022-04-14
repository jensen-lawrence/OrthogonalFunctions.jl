# ----------------------------------------------------------------------------------------------------------------------
# Chebyshev Polynomials of the First Kind
# ----------------------------------------------------------------------------------------------------------------------

"""
    chebyshev_T_f(n, x)

Evaluates the nth Chebyshev polymomial of the first kind at x ∈ ℝ.
"""
function chebyshev_T_f(n::Int, x::R) where {R<:Real}
    T = sum([genbinom(n, 2k) * (x^2 - 1)^k * x^(n - 2k) for k ∈ 0:floor(Int, n/2)])
    return T
end

"""
    chebyshev_T_w(x)

Evaluates the orthogonality weight function of the nth Chebyshev polynomial 
of the first kind at x ∈ ℝ.
"""
function chebyshev_T_w(x::R) where {R<:Real}
    return √(1/(1 - x^2))
end

"""
    ChebyshevT(n)

Representation of the nth Chebyshev polynomial of the first kind, Tₙ, where
n ∈ ℤ⁺.

Available methods:
- `ChebyshevT(n).n` returns n (the polynomial degree).
- `ChebyshevT(n).f` returns the polyomial function.
- `ChebyshevT(n).w` returns the orthogonality weight function.
- `ChebyshevT(n).I` returns the orthogonality interval.

`ChebyshevT(n)` is directly callable, so Tₙ(x) can be evaluated using 
`ChebyshevT(n)(x)`. Tₙ(x) can also be evaluated using`ChebyshevT(n).f(x)`.
Similarly, the orthogonality weight function w(x) can be evaluated using 
`ChebyshevT(n).w(x)`.

Aliases: `T`, `ChebyshevFirst`.
"""
struct ChebyshevT <: AbstractChebyshev 
    n::Int 
    f
    w 
    I 

    function ChebyshevT(n)
        n ≥ 0 || error("n must be non-negative.")
        f = x -> chebyshev_T_f(n, x)
        w = x -> chebyshev_T_w(x)
        I = (-1, 1)
        new(n, f, w, I)
    end
end

# Make ChebyshevT(n) callable 
function (T::ChebyshevT)(x::R) where {R<:Real}
    return T.f(x)
end

"""
    T(n)

Alias of `ChebyshevT(n)`.
"""
function T(n::Int)
    return ChebyshevT(n)
end

"""
    ChebyshevFirst(n)

Alias of `ChebyshevT(n)`.
"""
function ChebyshevFirst(n::Int)
    return ChebyshevT(n)
end

"""
    innerproduct(T1, T2)

Returns the inner product of `T1` and `T2` using the orthogonality relation 
for Chebyshev polynomials of the first kind.
"""
function innerproduct(T1::ChebyshevT, T2::ChebyshevT)
    n, m = T1.n, T2.n 
    if n == m == 0
        return π
    else
        return π/2 * δ(n, m)
    end
end

# ----------------------------------------------------------------------------------------------------------------------
# Chebyshev Polynomials of the Second Kind
# ----------------------------------------------------------------------------------------------------------------------

"""
    chebyshev_U_f(n, x)

Evaluates the nth Chebyshev polymomial of the second kind at x ∈ ℝ.
"""
function chebyshev_U_f(n::Int, x::R) where {R<:Real}
    T = sum([genbinom(n + 1, 2k + 1) * (x^2 - 1)^k * x^(n - 2k) for k ∈ 0:floor(Int, n/2)])
    return T
end

"""
    chebyshev_U_w(x)

Evaluates the orthogonality weight function of the nth Chebyshev polynomial 
of the second kind at x ∈ ℝ.
"""
function chebyshev_U_w(x::R) where {R<:Real}
    return √(1 - x^2)
end

"""
    ChebyshevU(n)

Representation of the nth Chebyshev polynomial of the second kind, Uₙ, where
n ∈ ℤ⁺.

Available methods:
- `ChebyshevU(n).n` returns n (the polynomial degree).
- `ChebyshevU(n).f` returns the polyomial function.
- `ChebyshevU(n).w` returns the orthogonality weight function.
- `ChebyshevU(n).I` returns the orthogonality interval.

`ChebyshevU(n)` is directly callable, so Uₙ(x) can be evaluated using 
`ChebyshevU(n)(x)`. Uₙ(x) can also be evaluated using`ChebyshevU(n).f(x)`.
Similarly, the orthogonality weight function w(x) can be evaluated using 
`ChebyshevU(n).w(x)`.

Aliases: `U`, `ChebyshevSecond`.
"""
struct ChebyshevU <: AbstractChebyshev 
    n::Int 
    f
    w 
    I 

    function ChebyshevU(n)
        n ≥ 0 || error("n must be non-negative.")
        f = x -> chebyshev_U_f(n, x)
        w = x -> chebyshev_U_w(x)
        I = (-1, 1)
        new(n, f, w, I)
    end
end

# Make ChebyshevU(n) callable 
function (U::ChebyshevU)(x::R) where {R<:Real}
    return U.f(x)
end

"""
    U(n)

Alias of `ChebyshevU(n)`.
"""
function U(n::Int)
    return ChebyshevU(n)
end

"""
    ChebyshevSecond(n)

Alias of `ChebyshevU(n)`.
"""
function ChebyshevSecond(n::Int)
    return ChebyshevU(n)
end

"""
    innerproduct(U1, U2)

Returns the inner product of `U1` and `U2` using the orthogonality relation 
for Chebyshev polynomials of the second kind.
"""
function innerproduct(U1::ChebyshevU, U2::ChebyshevU)
    n, m = U1.n, U2.n 
    return π/2 * δ(n, m)
end

# ----------------------------------------------------------------------------------------------------------------------
# Chebyshev Polynomials of the Third Kind
# ----------------------------------------------------------------------------------------------------------------------

"""
    chebyshev_V_f(n, x)

Evaluates the nth Chebyshev polymomial of the third kind at x ∈ ℝ.
"""
function chebyshev_V_f(n::Int, x::R) where {R<:Real}
    if n ≥ 1
        return ChebyshevU(n).f(x) - ChebyshevU(n - 1).f(x)
    else
        return 1
    end
end

"""
    chebyshev_V_w(x)

Evaluates the orthogonality weight function of the nth Chebyshev polynomial 
of the third kind at x ∈ ℝ.
"""
function chebyshev_V_w(x::R) where {R<:Real}
    return √((1 + x)/(1 - x))
end

"""
    ChebyshevV(n)

Representation of the nth Chebyshev polynomial of the third kind, Vₙ, where
n ∈ ℤ⁺.

Available methods:
- `ChebyshevV(n).n` returns n (the polynomial degree).
- `ChebyshevV(n).f` returns the polyomial function.
- `ChebyshevV(n).w` returns the orthogonality weight function.
- `ChebyshevV(n).I` returns the orthogonality interval.

`ChebyshevV(n)` is directly callable, so Vₙ(x) can be evaluated using 
`ChebyshevV(n)(x)`. Vₙ(x) can also be evaluated using`ChebyshevV(n).f(x)`.
Similarly, the orthogonality weight function w(x) can be evaluated using 
`ChebyshevV(n).w(x)`.

Aliases: `V`, `ChebyshevThird`.
"""
struct ChebyshevV <: AbstractChebyshev 
    n::Int 
    f
    w 
    I 

    function ChebyshevV(n)
        n ≥ 0 || error("n must be non-negative.")
        f = x -> chebyshev_V_f(n, x)
        w = x -> chebyshev_V_w(x)
        I = (-1, 1)
        new(n, f, w, I)
    end
end

# Make ChebyshevV(n) callable 
function (V::ChebyshevV)(x::R) where {R<:Real}
    return V.f(x)
end

"""
    V(n)

Alias of `ChebyshevV(n)`.
"""
function V(n::Int)
    return ChebyshevV(n)
end

"""
    ChebyshevThird(n)

Alias of `ChebyshevV(n)`.
"""
function ChebyshevThird(n::Int)
    return ChebyshevV(n)
end

"""
    innerproduct(V1, V2)

Returns the inner product of `V1` and `V2` using the orthogonality relation 
for Chebyshev polynomials of the third kind.
"""
function innerproduct(V1::ChebyshevV, V2::ChebyshevV)
    n, m = V1.n, V2.n 
    return W * δ(n, m)
end

# ----------------------------------------------------------------------------------------------------------------------
# Chebyshev Polynomials of the Fourth Kind
# ----------------------------------------------------------------------------------------------------------------------

"""
    chebyshev_W_f(n, x)

Evaluates the nth Chebyshev polymomial of the fourth kind at x ∈ ℝ.
"""
function chebyshev_W_f(n::Int, x::R) where {R<:Real}
    if n ≥ 1
        return ChebyshevU(n).f(x) + ChebyshevU(n - 1).f(x)
    else
        return 1
    end
end

"""
    chebyshev_W_w(x)

Evaluates the orthogonality weight function of the nth Chebyshev polynomial 
of the third kind at x ∈ ℝ.
"""
function chebyshev_W_w(x::R) where {R<:Real}
    return √((1 - x)/(1 + x))
end

"""
    ChebyshevW(n)

Representation of the nth Chebyshev polynomial of the fourth kind, Wₙ, where
n ∈ ℤ⁺.

Available methods:
- `ChebyshevW(n).n` returns n (the polynomial degree).
- `ChebyshevW(n).f` returns the polyomial function.
- `ChebyshevW(n).w` returns the orthogonality weight function.
- `ChebyshevW(n).I` returns the orthogonality interval.

`ChebyshevW(n)` is directly callable, so Wₙ(x) can be evaluated using 
`ChebyshevW(n)(x)`. Wₙ(x) can also be evaluated using`ChebyshevW(n).f(x)`.
Similarly, the orthogonality weight function w(x) can be evaluated using 
`ChebyshevW(n).w(x)`.

Aliases: `W`, `ChebyshevFourth`.
"""
struct ChebyshevW <: AbstractChebyshev 
    n::Int 
    f
    w 
    I 

    function ChebyshevW(n)
        n ≥ 0 || error("n must be non-negative.")
        f = x -> chebyshev_W_f(n, x)
        w = x -> chebyshev_W_w(x)
        I = (-1, 1)
        new(n, f, w, I)
    end
end

# Make ChebyshevW(n) callable 
function (W::ChebyshevW)(x::R) where {R<:Real}
    return W.f(x)
end

"""
    W(n)

Alias of `ChebyshevW(n)`.
"""
function W(n::Int)
    return ChebyshevW(n)
end

"""
    ChebyshevFourth(n)

Alias of `ChebyshevW(n)`.
"""
function ChebyshevFourth(n::Int)
    return ChebyshevW(n)
end

"""
    innerproduct(W1, W2)

Returns the inner product of `W1` and `W2` using the orthogonality relation 
for Chebyshev polynomials of the fourth kind.
"""
function innerproduct(W1::ChebyshevW, W2::ChebyshevW)
    n, m = W1.n, W2.n 
    return π * δ(n, m)
end

# ----------------------------------------------------------------------------------------------------------------------