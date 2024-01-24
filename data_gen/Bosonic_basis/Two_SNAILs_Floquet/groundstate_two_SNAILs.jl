include("../../../src/main_2_snails.jl") 


using DataFrames # this is like pandas
using CSV 
using ProgressBars


N = 15  # this is so dim(H_a⊗H_b) = NxN. i.e. if N=10 ⟹ dim(H_a⊗H_b) = 100.

# photon number operators. useful for filtering data later on.
n_op_a = a(N)'*a(N)
n_op_b = b(N)'*b(N)

##############
# Parameters #
##############

# ALL ENERGY/FREQUENCY (ħ=1) UNITS ARE IN kHz !!!!!!!!!!

# Undriven part parameters
ω_0_a = (2*π) * 5.97*1e6 # kHz
ω_0_b = (2*π) * 5.9888*1e6 # kHz
g_3_a, g_4_a = 3*(2*π) * 5e4, 4*(2*π) * 5e2 # kHz trust that is multiplied by 3,4.
g_3_b, g_4_b = 3*(2*π) * 5e3, 4*(2*π) * 5e1 # kHz
# g_3_a, g_4_a = 3*(2*π) * 5e3, 4*(2*π) * 5e1 # kHz trust that is multiplied by 3,4.
# g_3_b, g_4_b = 3*(2*π) * 5e3, 4*(2*π) * 5e1 # kHz
g_array = (2*π) * 1e4 .* Vector(range(0,1,length=100)) # kHz


# Drive parameters
Ω_a = (2*π) * 1e6 # kHz
Ω_b = (2*π) * 5e5 # kHz

ω_a = 2*ω_0_a # in kHz
ω_b = 2*ω_0_b # in kHz

# Derived parameters 
K_a = (10*g_3_a^2)/(3*ω_0_a) - 3*g_4_a/2
K_b = (10*g_3_b^2)/(3*ω_0_b) - 3*g_4_b/2

ϵ_2_a = g_3_a * Π(ω_0_a,Ω_a,ω_a)  
ϵ_2_b = g_3_b * Π(ω_0_b,Ω_b,ω_b)  


ω_a_n, ω_b_m = n_m_conmens_period(Int(round(ω_a/(2*π),sigdigits=6)),Int(round(ω_b/(2*π),sigdigits=6)))

# T = ω_a_n * (2*π)/(ω_a) # in kHz^-1
T = (2*π)/ω_b  # this is in case you want to hard code it, the period is very large eh.


################
# Floquet code #
################

# Define data form: ΔE_n | J | n_photons_a | n_photons_b
floquet_data_array = zeros( (N^2)*length(g_array), 4 )

# Propose a coherent superposition as initial groundstate of each well
α_A = √(Complex(ϵ_2_a/K_a))
ψ_0_a = coherent_state(N,α_A) + coherent_state(N,-α_A)
ψ_0_a = ψ_0_a/norm(ψ_0_a)

α_B = √(Complex(ϵ_2_b/K_b))
ψ_0_b = coherent_state(N,α_B) + coherent_state(N,-α_B)
ψ_0_b = ψ_0_b/norm(ψ_0_b)

V_0 = kron(ψ_0_a,ψ_0_b)

counter = 1
pbar = ProgressBar(total=length(g_array))
for (j,g) in enumerate(g_array)
    # find closest groundstate and eigenergies and update
    ϵ_n, ϕ_n = qen_qmodes_2_snails(N,T, g, ω_0_a, ω_0_b, g_3_a, g_4_a, g_3_b, g_4_b, Ω_a, ω_a, Ω_b, ω_b)

    overlaps = zeros(N^2)
    for n=1:N^2
        overlaps[n] = real( (V_0'*ϕ_n[:,n])*(ϕ_n[:,n]'*V_0) )
    end
    ϵ_0 = ϵ_n[argmax(overlaps)]
    V_0 = ϕ_n[:,argmax(overlaps)]

    # Δϵ_n = sort( mod.( ϵ_n .- ϵ_0, ω_2/2) )
    Δϵ_n =  mod.( ϵ_n .- ϵ_0, 2*π/T) 

    #store parameters in array
    for n=1:N^2
        n_photons_a = real( ϕ_n[:,n]'*n_op_a*ϕ_n[:,n] )
        n_photons_b = real( ϕ_n[:,n]'*n_op_b*ϕ_n[:,n] )

        floquet_data_array[counter,:] .= Δϵ_n[n], g, n_photons_a, n_photons_b
        counter += 1
    end
    update(pbar)
end

# Put data in convenient DataFrame object and save it
df_floquet_2_snails = DataFrame(floquet_data_array, ["ΔE_n","g","n_photons_a","n_photons_b"]) 
# CSV.write("data/floquet_crossings_comparison.csv", df_floquet)

# Combine both sims into a single dataframe
# df_comparison = DataFrame(vcat(floquet_data_array,H_eff_data_array), ["ΔE_n","ϵ_1","ϵ_2","n_photons","Floquet"]) 

# CSV.write("data/crossings_comparison.csv", df_comparison)

display(norm(kron(ψ_0_a,ψ_0_b)'*V_0))