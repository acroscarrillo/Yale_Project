using QuantumOptics # to plot wigner functions
using Roots # to find minimas 
using Flux # to find gradients for minimas 
using Test # to throw warning

# utility code goes here
function coherent_state(N,α)
    A = a(N)
    vac = zeros(N)
    vac[1] = 1
    coh_state =  exp(-α*A+(α*A)')*vac
    return coh_state/norm(coh_state)
end

function coherent_state(N,x,p)
    α =  (x + im*p)/√(2) 
    return coherent_state(N,α)
end

function IPR(vec)
    return sum( [norm(vec[n])^4 for n=1:length(vec)] )
end

function which_basis_state(v,basis_matrix)
    basis_overlap = basis_matrix' * v # check multiplication is correct
    max_value, index = findmax(basis_overlap)
    return max_value, index
end

function avg_spacing_ratio(Δnn_array)
    N = length(Δnn_array) - 1

    r_tidle = 0
    for i=1:N
        r = Δnn_array[i+1]/Δnn_array[i]    
        r_tidle += r/N
    end

    return r_tidle
end

function avg_spacing_ratio_min(Δnn_array)
    N = length(Δnn_array) - 1

    r_tidle = 0
    for i=1:N
        r = Δnn_array[i+1]/Δnn_array[i]
        if r < 1/r
            r_tidle += r/N
        else
            r_tidle += 1/(r*N)
        end
    end
    return r_tidle
end

function rnm_avg_spacing_ratio_min(Δnn_array)
    r_coe = 0.53
    r_p = 0.39
    numerator = avg_spacing_ratio_min(Δnn_array) - r_p
    denominator = r_coe - r_p
    return numerator/denominator
end

function rand_sym(N)
    rand_mat = randn(N,N) 
    return (rand_mat + rand_mat')
end

function rand_herm(N)
    rand_mat = randn(N,N) + im*randn(N,N) 
    return (rand_mat + rand_mat')
end

function dblwell_minimas(K,ϵ_1,ϵ_2,xlim=-20) # 20 should be fine?
    H_temp(x) = H_cl(x,0,K,ϵ_1,ϵ_2)
    H_temp_first_derivative(x) = gradient(H_temp, x)[1]
    H_temp_second_derivative(x) = gradient(H_temp_first_derivative, x)[1]
    extremas =  find_zeros(H_temp_first_derivative,-xlim,xlim)
    minimas = []
    for x in extremas
        if H_temp_second_derivative(x) > 0
            push!(minimas,x)
        end
    end
    return sort(minimas,rev=true) #so that deep well is at index [1]
end

function wigner_func(ψ,xlim,plim,meshstep=0.1)
    N = length(ψ)
    basis = FockBasis(N-1)
    ψ_ket = Ket(basis,ψ)
    wigner_f = wigner(ψ_ket,-xlim:meshstep:xlim,-plim:meshstep:plim)
    return wigner_f' #transpose so that x is in x-axis
end

# Spin  stuff
function decimal_to_bit(L,n)
    if n > (2^L-1)
        max_fock = 2^L-1
        throw(ArgumentError("Number state |$n> cannot be represented in a spin chain of length L=$L. You may at most represent |$max_fock>."))
    end
    bitstring_rep = bitstring(Int32(n))
    bit_array = zeros(L)
    for (i,character) in enumerate( bitstring_rep[end-L+1:end] )
        bit_array[i] = parse(Int64,character)
    end
    return bit_array
end

function fock_basis_to_spin_basis(L,n_ocupation)
    spin_basis = zeros(2^L)
    spin_basisn_ocupation
    return spin_basis
end

# function fock_op_to_spin_op(Op)
#     L = size(Op)[1]
#     if L > 18
#         Warn("Warning! spin chain exceeds L=18 with length $L.")
#     end

#     spin_op = zeros(2^L,2^L)
#     for i=1:L
#         for j=1:L

        

    
# end

function anharmonicity(H)
    lamb, _ = eigen(H)
    return lamb[3]+lamb[1]-2*lamb[2]
end
