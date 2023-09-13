using LinearAlgebra # to diagonalise stuff
using DifferentialEquations


"""
    a(N::Int)

Return the bosonic annahilation operator truncated at dimension `N`.
For the corresponding bosonic creation operator simply make use of the native adjoint operator `a(N)'`.

# Examples
```julia-repl
julia> a(3)
3×3 Matrix{Float64}:
 0.0  1.0  0.0
 0.0  0.0  1.41421
 0.0  0.0  0.0

julia> a(3)'
3×3 adjoint(::Matrix{Float64}) with eltype Float64:
 0.0  0.0      0.0
 1.0  0.0      0.0
 0.0  1.41421  0.0
```
"""
function a(N::Int)
    a = zeros(Float64, (N, N))
    for n in 1:N-1
        a[n,n+1] = sqrt(n)
    end
    return a
end


"""
    H(N,ω_0,g_n,Ω_d,ω_d,t)

Return the experimental Kerr cat Hamiltonian as given in, for instance, in Eq. (1) of https://arxiv.org/pdf/2209.03934.pdf. Non-liniarities `g_n` should be passed as vector whose componets are assigned as `g_3 = g_n[1]`, `g_4 = g_n[2]`, etc. 

# Examples
```julia-repl
julia> H(3,1,[1,1],1,2,1)
3×3 Matrix{ComplexF64}:
     3.0+0.0im           3.0+0.416147im  4.24264+0.0im
     3.0-0.416147im     10.0+0.0im       4.24264+0.588521im
 4.24264+0.0im       4.24264-0.588521im      8.0+0.0im

```
"""
function H(N,ω_0,g_n,Ω_d,ω_d,t)
    A = a(N) # annahilation op. up to dim N
    exp_order = length(g_n)
    expansion = [(A' + A)^n for n=3:(2+exp_order)]
    return  ω_0*A'*A + g_n'*expansion - im*Ω_d*cos(ω_d*t)*(A - A')
end


"""
    H_0(N,ω_0,g_n)

Return the time-independent part of the experimental Kerr cat Hamiltonian as given in, for instance, Eq. (1) of https://arxiv.org/pdf/2209.03934.pdf. Non-liniarities `g_n` should be passed as vector whose componets are assigned as `g_3 = g_n[1]`, `g_4 = g_n[2]`, etc. This function is equivalent to `H(N,ω_0,g_n,0,0,0)`. Note the full Hamiltonian would be `H = H_0 + H_d*cos(ω_d*t)`.

# Examples
```julia-repl
julia> H_0(3,1,[1,1])
3×3 Matrix{Float64}:
 0.75     1.0      1.06066
 1.0      3.25     1.41421
 1.06066  1.41421  3.5

```
"""
function H_0(N,ω_0,g_n)
    A = a(N) # annahilation op. up to dim N
    exp_order = length(g_n)
    expansion = [(A' + A)^n for n=3:(2+exp_order)]
    return  ω_0*A'*A + g_n'*expansion
end



"""
    H_d(N)

Return the operator component of the time-dependent part of the experimental Kerr cat Hamiltonian as given in, for instance, Eq. (1) of https://arxiv.org/pdf/2209.03934.pdf. This is 

`- 2*im*(a - a')`

Note this term does NOT depend on time, the full Hamiltonian would be `H = H_0 + H_d*Ω_d*cos(ω_d*t)` if it only has one drive.

# Examples
```julia-repl
julia> H_d(3)
3×3 Matrix{ComplexF64}:
  0.0-0.0im   0.0-2.0im      0.0-0.0im
 -0.0+2.0im   0.0-0.0im      0.0-2.82843im
  0.0-0.0im  -0.0+2.82843im  0.0-0.0im

```
"""
function H_d(N)
    A = a(N) # annahilation op. up to dim N
    return  - 2*im*(A - A')
end

"""
    H_eff(N,Δ,K,ϵ_1,ϵ_2)

Return the effective Hamiltonian of the experimental Kerr cat Hamiltonian as explained in, for instance, Eq. (8) of https://arxiv.org/pdf/2210.07255.pdf, where however, a self-energy and a secondary driving term are missing (perhaps there's a better ref of this?).This is 

`Δ*a'*a - K*(a'^2)*(a^2) + ϵ_1*(a + a') + ϵ_2*(a^2 + a'^2)`

# Examples
```julia-repl
julia> H_eff(3,1,2,3,4)
3×3 Matrix{Float64}:
 -4.0      3.0      5.65685
  3.0      1.0      4.24264
  5.65685  4.24264  2.0

```
"""
function H_eff(N,Δ,K,ϵ_1,ϵ_2)
    A = a(N) # annahilation op. up to dim N
    return Δ*A'*A - K*(A'^2)*(A^2) + ϵ_1*(A + A') + ϵ_2*(A^2 + A'^2)
end

"""
    f!(u, p, t)

In-place function used to initiate the function `ODEProblem` from `DifferentialEquations.jl` that defines the propagator ODE problem

`d/dt U(t) = -i H(t) U(t),`

where `f! =  -i H(t) U(t)`. This will be passed to the function `solver` from `DifferentialEquations.jl` to obtain `U(t)` in functions like `quasienergies` and `qen_qmodes` to get the Floquet quasienergies and quasimodes. Please consult `DifferentialEquations.jl` documentation, in particular https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction. In our present case we have 

`-i*H(t) = -i*( H_0 + H_d*( Ω_1*cos(ω_1*t) + Ω_2*cos(ω_2*t) ) )`

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
    p[7] .= (-im) .* (p[1] .+ p[2] .* (p[3] .* cos(p[4]*t) .+ p[5] .* cos(p[6]*t)) )
    mul!(du, p[7], u) 
end

"""
    quasienergies(N, ω_0, g_n, Ω_d, ω_d)

Return the Floquet quasienergies of H(N,ω_0,g_n,Ω_d,ω_d,t). Needs better documentation...

# Warning!!
The period of the system is hard coded to be `2*pi/ω_2`, but this is not true in general! This is so because in the experiment, ω_0 ≈ ω_1 ≈ ω_2/2, so ω_2 will have the largest period (these guys are commensurable).
"""
function quasienergies(N, ω_0, g_n, Ω_1, ω_1, Ω_2, ω_2)
    p =  H_0(N,ω_0,g_n), H_d(N), Ω_1, ω_1, Ω_2, ω_2, ComplexF64.(zeros(N,N))
    T = 2*pi/ω_2 #THIS IS NOT GENERAL, see documentation above.
    tspan = (0.0, 2*T) # because rodrigo said so (need to understand this)
    u_0 = ComplexF64.(Matrix(I(N)))
    prob = ODEProblem(f!, u_0, tspan, p)
    sol = solve(prob) 
    U_2T = sol.u[end]
    η_n, _ = eigen(U_2T)
    ϵ_n = im*log.(η_n)/T
    return mod.(real.(ϵ_n), ω_2/2)  # as they are mod(ω_d/2), cause we take T to be 2T.
end

"""
    qen_qmodes(N, ω_0, g_n, Ω_1, ω_1, Ω_2, ω_2)

Return the Floquet quasienergies and quasimodes of H(N,ω_0,g_n,Ω_d,ω_d,t).Needs better documentation...

# Warning!!
The period of the system is hard coded to be `2*pi/ω_2`, but this is not true in general! This is so because in the experiment, ω_0 ≈ ω_1 ≈ ω_2/2, so ω_2 will have the largest period (these guys are commensurable).
"""
function qen_qmodes(N, ω_0, g_n, Ω_1, ω_1, Ω_2, ω_2)
    p =  H_0(N,ω_0,g_n), H_d(N), Ω_1, ω_1, Ω_2, ω_2, ComplexF64.(zeros(N,N))
    T = 2*pi/ω_2 #THIS IS NOT GENERAL, see documentation above.
    tspan = (0.0, 2*T) # because rodrigo said so (need to understand this)
    u_0 = ComplexF64.(Matrix(I(N)))
    prob = ODEProblem(f!, u_0, tspan, p)
    sol = solve(prob) 
    U_2T = sol.u[end]
    η_n, ϕ_n = eigen(U_2T)
    ϵ_n = im*log.(η_n)/T
    return mod.(real.(ϵ_n), ω_2/2), ϕ_n  # as they are mod(ω_d/2)
end