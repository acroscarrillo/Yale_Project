using LinearAlgebra 
using DifferentialEquations

"""
    Q(N::Int)

Return the charge operator `Q` truncated at dimension `N`. Please note that `N` must be a an odd number greater than 1.

# Examples
```julia-repl
julia> Q(5)
5×5 Matrix{Int64}:
 -2   0  0  0  0
  0  -1  0  0  0
  0   0  0  0  0
  0   0  0  1  0
  0   0  0  0  2
```
"""
function Q(N::Int)
    if iseven(N)||N<3
        throw(ArgumentError("Hilbert space dimension, N, must be an odd number greater than 1 but N=$N was given."))
    end
    return Matrix(diagm(-N÷2:N÷2))
end


"""
    exp_iϕ(N::Int)

Return the flux operator `exp_iϕ` truncated at dimension `N`. We chose to enforce `N` to be an odd number greater than 1 just to match this same requirement on `Q`. Note that `exp_iϕ` is real unitary.

# Examples
```julia-repl
julia> exp_iϕ(5)
5×5 Matrix{Float64}:
 0.0  1.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0
 0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  1.0
 1.0  0.0  0.0  0.0  0.0
```
"""
function exp_iϕ(N::Int)
    if iseven(N)||N<3
        throw(ArgumentError("Hilbert space dimension, N, must be an odd number greater than 1 but N=$N was given."))
    end
    temp_mat = diagm(1=>ones(N-1))
    temp_mat[end,1] = 1
    return temp_mat
end


"""
    cos_ϕ(N,ϕ_ext,m)

Returns the operator 

`` cos( (ϕ_ext-ϕ)/m ) ``

truncated at dimension `N`. This is done by decomposing the cosine in complex exponentials so that we can call `exp_iϕ`. That is,

`` cos( (ϕ_ext-ϕ)/m ) = 1/2 * ( exp(iϕ_ext/m)exp(-iϕ/m)  + exp(-iϕ_ext/m)exp(iϕ/m)) ``,

although `exp(-iϕ/m)` is simply the transpose of `exp(iϕ/m)` as it's unitary.

# Arguments
    - `N`: Hilbert space dimension, an integer number.
    - `ϕ_ext`: External flux, a real number (not an operator!).
    - `m`: An integer number.

# Examples
```julia-repl
julia> round.( cos_ϕ(3,1,3), sigdigits=2 )
3×3 Matrix{ComplexF64}:
   0.8+0.0im   0.074+0.12im  0.074-0.12im
 0.074-0.12im    0.8+0.0im   0.074+0.12im
 0.074+0.12im  0.074-0.12im    0.8+0.0im
```
"""
function cos_ϕ(N,ϕ_ext,m)
    exp_tmp = exp_iϕ(N)^(1/m)
    return 0.5*( exp(im*ϕ_ext/m)*exp_tmp + exp(-im*ϕ_ext/m)*exp_tmp' )
end


"""
    H_joseph(N,α,E_J,ϕ_ext,m)

Returns the operator Josephson Hamiltonian

``α*E_J*cos(ϕ) + m*E_J*cos( (ϕ_ext-ϕ)/m )``

# Examples
```julia-repl
julia> round.( H_joseph(3,1,5,7,3), sigdigits=2 )
3×3 Matrix{ComplexF64}:
 -8.7+0.0im   1.7+4.0im   1.7-4.0im
  1.7-4.0im  -8.7+0.0im   1.7+4.0im
  1.7+4.0im   1.7-4.0im  -8.7+0.0im
```
"""
function H_joseph(N,α,E_J,ϕ_ext,m)
    return α*E_J*cos_ϕ(N,0,1) + m*E_J*cos_ϕ(N,ϕ_ext,m)
end