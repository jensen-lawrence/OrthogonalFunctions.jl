"""
    fac(x)

Generalized factorial function. Returns:
- `factorial(x)` if 0 ≤ x ≤ 20 and x ∈ ℤ⁺
- `gamma(x + 1)` otherwise
"""
function fac(x::T) where {T <: Real}
    if (typeof(x) <: Integer) && (0 ≤ x ≤ 20)
        return factorial(x)
    else
        return gamma(x + 1)
    end
end

"""
    binom(x, y)

Generalized binomial coefficient. Returns:
- `binomial(x, y)` if x, y ∈ ℤ⁺ and 0 ≤ x ≤ 66 
- `gamma(x + 1)/(gamma(y + 1) * gamma(x - y + 1))` otherwise
"""
function binom(x::T1, y::T2) where {T1, T2 <: Real}
    if (typeof(x) <: Integer) && (typeof(y) <: Integer) && (0 ≤ x ≤ 66)
        return binomial(x, y)
    else
        return gamma(x + 1)/(gamma(y + 1) * gamma(x - y + 1))
    end
end

"""
    δ(n, m)

Kronecker delta function.
"""
function δ(n::Integer, m::Integer)
    if n == m 
        return 1 
    else 
        return 0
    end
end