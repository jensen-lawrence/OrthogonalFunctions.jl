"""
    genfac(x)

Generalized factorial function. Cases:
- returns `factorial(x)` if 0 ≤ x ≤ 20 and x ∈ ℤ⁺
- returns `gamma(x + 1)` otherwise
"""
function genfac(x::T) where {T<:Real}
    if (typeof(x) <: Int) && (0 ≤ x ≤ 20)
        return factorial(x)
    else
        return gamma(x + 1)
    end
end

"""
    genbinom(x, y)

Generalized binomial coefficient. Cases:
- returns `binomial(x, y)` if x, y ∈ ℤ⁺ and 0 ≤ x ≤ 66
- returns `gamma(x + 1)/(gamma(y + 1) * gamma(x - y + 1))` otherwise 
"""
function genbinom(x::T1, y::T2) where {T1,T2<:Real}
    if (typeof(x) <: Int) && (typeof(y) <: Int) && (0 ≤ x ≤ 66)
        return binomial(x, y)
    else
        return gamma(x + 1)/(gamma(y + 1)*gamma(x - y + 1))
    end
end

"""
    δ(n, m)

Kronecker delta function.
"""
function δ(n::Int, m::Int)
    if n == m
        return 1
    else
        return 0
    end
end