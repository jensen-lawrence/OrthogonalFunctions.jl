module Orthogonomial

# Imports
using Roots: find_zeros
using SpecialFunctions: gamma 
using LaTeXStrings: latexstring
using Latexify: latexify
using SymPy

# Files
include("types.jl")
include("utils.jl")
include("polynomials/polynomials.jl")

end