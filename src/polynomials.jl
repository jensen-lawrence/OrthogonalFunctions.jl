"""
    degree(F)

Returns the degree of `F`.
"""
degree(F::OrthogonalPolynomial) = F.n 

"""
    polynomial(F)

Returns the polynomial function of `F`.
"""
polymomial(F::OrthogonalPolynomial) = F.f 

"""
    weight(F)

Returns the orthogonality weight function of `F`.
"""
weight(F::OrthogonalPolynomial) = F.w 

"""
    interval(F)

Returns the orthogonality interval of `F`.
"""
interval(F::OrthogonalPolynomial) = F.I 