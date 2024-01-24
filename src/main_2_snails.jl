using LinearAlgebra # to diagonalise stuff
using DifferentialEquations
using Primes
using Roots # to find minimas 
using Flux # to find gradients for minimas 
using QuantumOptics # to plot wigner functions


# define general anihilation
function c(N::Int)
    temp = zeros(Float64, (N, N))
    for n in 1:N-1
        temp[n,n+1] = sqrt(n)
    end
    return temp
end

function a(N::Int)
    temp = c(N)
    return kron(temp,I(N))
end

function b(N::Int)
    temp = c(N)
    return kron(I(N),temp)
end


function H_2_snails(N,g,ω_0_a,ω_0_b,g_3_a,g_4_a,g_3_b,g_4_b)
    A, B = a(N), b(N) # annahilation op. up to dim N
    term_self =  ω_0_a*A'*A + ω_0_b*B'*B
    term_3 = (g_3_a/3)*(A+A')^3 + (g_3_b/3)*(B+B')^3
    term_4 =  (g_4_a/4)*(A+A')^4 + (g_4_b/4)*(B+B')^4
    term_int = g*(A'*B + A*B') # A and B commute as they are in different tensor spaces
    return term_self + term_3 + term_4 + term_int
end

function H_d_a(N)
    A = a(N) # annahilation op. up to dim N
    return  - 1*im*(A - A')
end


function H_d_b(N)
    B = b(N) # annahilation op. up to dim N
    return  - 1*im*(B - B')
end




function n_m_conmens_period(int_1::Int,int_2::Int)
    if int_1<int_2
        prime_array = reverse( primes(0,int_1) )
        for p in prime_array
            while  int_1%p == 0 && int_2%p == 0 
                int_1 = int_1/p
                int_2 = int_2/p
            end
        end
    else
        prime_array = reverse( primes(0,int_2) )
        for p in prime_array
            while  int_1%p == 0 && int_2%p == 0 
                int_1 = int_1/p
                int_2 = int_2/p
            end
        end
    end
    return int_1, int_2
end

function f!(du,u, p, t)
    # -i*H(t) = -i*( H + H_d_a*Ω_a*cos(ω_a*t) + H_d_b*Ω_b*cos(ω_b*t) ) )
    p[8] .= (+im).*( p[1] .+ p[2] .* ( p[3] .* cos(p[4]*t-π/2) ) .+ p[5] .* ( p[6] .* cos(p[7]*t-π/2) ) )
    mul!(du, p[8], u) 
end

function qen_qmodes_2_snails(N,T, g, ω_0_a, ω_0_b, g_3_a, g_4_a, g_3_b, g_4_b, Ω_a, ω_a, Ω_b, ω_b)
    p = H_2_snails(N,g,ω_0_a,ω_0_b,g_3_a,g_4_a,g_3_b,g_4_b), H_d_a(N), Ω_a, ω_a, H_d_b(N), Ω_b, ω_b, ComplexF64.(zeros(N^2,N^2))
    tspan = (0.0, T) 
    u_0 = ComplexF64.(Matrix(I(N^2)))
    prob = ODEProblem{true,SciMLBase.FullSpecialize}(f!, u_0, tspan, p)
    sol = solve(prob,saveat=[T]) 
    U_T = sol.u[end]
    η_n, ϕ_n = eigen(U_T)
    ϵ_n = im*log.(η_n)/T
    return mod.(real.(ϵ_n), 2*π/T), ϕ_n  # as they are mod(ω_d/2)
end



# utility code goes here
function Π(ω_0,Ω_2,ω_2)
    return (Ω_2*ω_2)/(ω_2^2-ω_0^2)
end

function coherent_state(N,α)
    A = c(N)
    vac = zeros(N)
    vac[1] = 1
    coh_state =  exp(-α*A+(α*A)')*vac
    return coh_state/norm(coh_state)
end
function coherent_state(N,x,p)
    α =  (x + im*p)/√(2) 
    return coherent_state(N,α)
end

function wigner_func(ψ,xlim,plim,meshstep=0.1)
    N = length(ψ)
    basis = FockBasis(N-1)
    ψ_ket = Ket(basis,ψ)
    wigner_f = wigner(ψ_ket,-xlim:meshstep:xlim,-plim:meshstep:plim)
    return wigner_f' #transpose so that x is in x-axis
end

function partial_trace(ψ)
    N = length(ψ)
    N_A = Int(√(N))
    projectors = []
    for n=1:N_A
        fock_state = zeros(N_A)
        fock_state[n] = 1
        push!( projectors, kron(I(N_A),fock_state) )
    end

    ρ = ψ*ψ'
    ρ_sub = zeros(N_A,N_A)
    for proj in projectors
        ρ_sub += proj' * ρ * proj
    end
    return ρ_sub
end

function entropy(ρ)
    λ_array, _ = eigen(ρ)
    temp = 0
    for λ in λ_array
        λ = real(λ)
        if λ > 0
            temp +=  -λ*log2(λ)
        end
    end
    return temp
end

function entanglement_entropy(ψ)
    ρ_A = partial_trace(ψ)
    return entropy(ρ_A)
end