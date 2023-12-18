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
    H_joseph(N,E_J,α,ϕ_ext,m)

Returns the Josephson Hamiltonian

`` α*E_J*cos(ϕ) + m*E_J*cos( (ϕ_ext-ϕ)/m ) ``

# Examples
```julia-repl
julia> round.( H_joseph(3,1,5,7,3), sigdigits=2 )
3×3 Matrix{ComplexF64}:
 -1.7+0.0im    2.3+0.81im   2.3-0.81im
  2.3-0.81im  -1.7+0.0im    2.3+0.81im
  2.3+0.81im   2.3-0.81im  -1.7+0.0im
```
"""
function H_joseph(N,E_J,α,ϕ_ext,m)
    return α*E_J*cos_ϕ(N,0,1) + m*E_J*cos_ϕ(N,ϕ_ext,m)
end


"""
    H_joseph(N,α,E_J,ϕ_ext,m)

Returns the Snail Hamiltonian

`` 4*E_C*Q(N)^2 - α*E_J*cos(ϕ) - m*E_J*cos( (ϕ_ext-ϕ)/m ) ``

# Minor performance comment
It constructs `Q(N)` twice instead of one. It shouldnt be a problem as `Q(N)` is diagonal but oh well.

# Examples
```julia-repl
julia> round.( H_snail(3,1,2,5,7,3), sigdigits=2 )
3×3 Matrix{ComplexF64}:
 13.0-0.0im  -4.2-4.0im  -4.2+4.0im
 -4.2+4.0im   8.7-0.0im  -4.2-4.0im
 -4.2-4.0im  -4.2+4.0im  13.0-0.0im
```
"""
function H_snail(N,E_C,E_J,α,ϕ_ext,m)
    return 4*E_C*Q(N)^2 - H_joseph(N,α,E_J,ϕ_ext,m) 
end


"""
    V_d(N)

Returns the operatorial part of the drive. In this case is the operator 

`` √(2)*Q ``
"""
function V_d(N)
    return √(2)*Q(N)
end

"""
    f!(u, p, t)

In-place function used to initiate the function `ODEProblem` from `DifferentialEquations.jl` that defines the propagator ODE problem

`d/dt U(t) = -i H(t) U(t),`

where `f! =  -i H(t) U(t)`. This will be passed to the function `solver` from `DifferentialEquations.jl` to obtain `U(t)` in functions like `quasienergies` and `qen_qmodes` to get the Floquet quasienergies and quasimodes. Please consult `DifferentialEquations.jl` documentation, in particular https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction. In our present case we have 

`H(t) =  H_0 + H_d*( Ω_1*cos(ω_1*t) + Ω_2*cos(ω_2*t) ) `

note there is no need for a minus sign anymore as we have found internal peace.
...
# Arguments
- `du`: derivative of propagator, i.e. `-i H(t) U(t)`.
- `u`: propagator.
- `p`: parameters. p =  H_0, H_d, Ω_1, ω_1, Ω_2, ω_2, du
- `t`: time.
...

# Performance
Several performance comments are in order. Function `f!` is called many many times by the solver which means that the faster `f!` runs the faster the solver will be. This is why `H_0`, `H_d`, `ω_1`,... etc. are passed as parameters so not to re-compute them again uncessarily.  Finally, to avoid uncessary memory re-allocations, we also pass `du` as a parameter forcing the solver store `du` in the same memory space which is much more efficient (this is achieved with the broadcasting assigment). Another comment is in order: there doesnt seem to be any performance benefit (nor in time nor in memory allocations) between creating two functions, say `f1!` and `f2!`, where `f1!` is speciallised to one single drive and `f2!` with two drives. This is probably because the operator `H_d` is shared as in the time dependent part is `~ H_d*(Ω_1*cos(ω_1) + Ω_2*cos(ω_2))`. So the only extra added complexity lies in computing two cosines and adding two numbers which is neglegible. 
"""
function f!(du,u, p, t)
    # -i*H(t) = -i*( H_0 + H_d*( Ω_1*cos(ω_1*t) + Ω_2*cos(ω_2*t) ) )
    p[7] .= (-im).*(p[1] .+ p[2] .* (p[3] .* cos(p[4]*t-π/2) .+ p[5] .* cos(p[6]*t-π/2)) )
    mul!(du, p[7], u) 
end

"""
    qen_qmodes(N, E_C, E_J, α, ϕ_ext, m, Ω_1, ω_1, Ω_2, ω_2)

Return the Floquet quasienergies and quasimodes of the driven Snail. This function calls `H_snail(N,E_C,E_J,α,ϕ_ext,m)` and `V_d(N)` as documented in this file. If nothing has changed since this documentation was written, these two functions are given by 

`H_snail = 4*E_C*Q(N)^2 - α*E_J*cos(ϕ) - m*E_J*cos( (ϕ_ext-ϕ)/m )`,

and

`V_d = √(2)*Q`,

which multiplies two `Ω_1*cos(ω_1 t)+Ω_1*cos(ω_1 t)` terms one for each drive. If in doubt, please check the source code of this functions and in particular of `f!` which is of fundamental importance.

# Warning!!
The period of the system is hard coded to be `2*pi/ω_2`, but this is not true in general! This is so because in the experiment, ω_0 ≈ ω_1 ≈ ω_2/2, so ω_2 will have the largest period (these guys are commensurable).
"""
function qen_qmodes(N, E_C, E_J, α, ϕ_ext, m, Ω_1, ω_1, Ω_2, ω_2)
    p =  H_snail(N,E_C,E_J,α,ϕ_ext,m), V_d(N), Ω_1, ω_1, Ω_2, ω_2, ComplexF64.(zeros(N,N))
    T_d2 = 2*pi/ω_2 #THIS IS NOT GENERAL, see documentation above.
    tspan = (0.0, 2*T_d2) 
    u_0 = ComplexF64.(Matrix(I(N)))
    prob = ODEProblem(f!, u_0, tspan, p)
    sol = solve(prob) 
    U_2T = sol.u[end]
    η_n, ϕ_n = eigen(U_2T)
    ϵ_n = im*log.(η_n)/(2*T_d2)
    return mod.(real.(ϵ_n), ω_2/2), ϕ_n  # as they are mod(ω_d/2)
end