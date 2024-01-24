include("../../../src/main_2_snails_dress.jl") 


using DataFrames # this is like pandas
using CSV 
using ProgressBars


N = 10  # this is so dim(H_a⊗H_b) = NxN. i.e. if N=10 ⟹ dim(H_a⊗H_b) = 100.

# photon number operators. useful for filtering data later on.
n_op_a = a(N)'*a(N)
n_op_b = b(N)'*b(N)

#######################
# OG basis Parameters #
#######################

# ALL ENERGY/FREQUENCY (ħ=1) UNITS ARE IN kHz !!!!!!!!!!

# Undriven part parameters
# ω_0_a_og = (2*π) * 5.97*1e6 # kHz
# ω_0_b_og = (2*π) * 5.9888*1e6 # kHz
ω_0_a_og = (2*π) * 4*1e6 # kHz
ω_0_b_og = (2*π) * 4.2*1e6 # kHz
g_3_a_og, g_4_a_og = 3*(2*π) * 5e4, 4*(2*π) * 5e2 # kHz trust that is multiplied by 3,4.
g_3_b_og, g_4_b_og = 3*(2*π) * 5e4, 4*(2*π) * 5e1 # kHz
# g_3_a, g_4_a = 3*(2*π) * 5e3, 4*(2*π) * 5e1 # kHz trust that is multiplied by 3,4.
# g_3_b, g_4_b = 3*(2*π) * 5e3, 4*(2*π) * 5e1 # kHz

g_array = (2*π) * 1e4 .* Vector(range(0,1,length=100)) # kHz
g_4_array = (2*π) * 4 * Vector(range(5e2,2e4,length=100))


# Derived parameters 
K_a_og = (10*g_3_a_og^2)/(3*ω_0_a_og) - 3*g_4_a_og/2
K_b_og = (10*g_3_b_og^2)/(3*ω_0_b_og) - 3*g_4_b_og/2

function Ω_new(ω_0,ω_d,g_3,g_4)
    num = 4 * ( 10*g_3^2/(3*ω_0) -3*g_4/2 ) * (ω_d^2-ω_0^2)
    denom = g_3*ω_d
    return num/denom
end


##############################################
# maping OG basis parameters to dress basis  #
##############################################
function new_params(N,g_og,ω_0_a_og,ω_0_b_og,g_3_a_og,g_4_a_og,g_3_b_og,g_4_b_og)
    Δ = Δ_func(ω_0_a_og,ω_0_b_og)
    θ = θ_func(g_og,ω_0_a_og,ω_0_b_og)
    #ω_0_a = 0.5 * ( ω_0_a_og + ω_0_b_og - √( Δ^2 + 4*g_og^2) )
    #ω_0_b =  0.5 * ( ω_0_a_og + ω_0_b_og + √( Δ^2 + 4*g_og^2) )

    ω_0_a = ω_0_a_og # change from Max
    ω_0_b = ω_0_b_og # change from Max

    g_3_a, g_3_b = g_3_a_og, g_3_b_og
    g_4_a, g_4_b = g_4_a_og, g_4_b_og

    ω_p = 0
    # ω_p = 2*ω_0_a
    ω_d = 2*ω_0_b

    # Ω_p = Ω_new(ω_0_a,ω_p,g_3_a,g_4_a)
    Ω_p = 0
    Ω_d = Ω_new(ω_0_b,ω_d,g_3_b,g_4_b)

    # display("ω_p=$ω_p")
    # display("ω_d=$ω_d")


    ω_p_n, ω_d_m = n_m_conmens_period(Int(round(ω_p/(2*π),sigdigits=6)),Int(round(ω_d/(2*π),sigdigits=6)))

    ϵ_2_a = g_3_a * Π(ω_0_a,Ω_p,ω_p)  
    ϵ_2_b = g_3_b * Π(ω_0_b,Ω_d,ω_d)  

    # T = ω_p_n * (2*π)/(ω_p) # in kHz^-1
    T = (2*π)/(ω_d)

    return N, T, θ, ω_0_a, ω_0_b, g_3_a, g_4_a, g_3_b, g_4_b, Ω_p, ω_p, Ω_d, ω_d, ϵ_2_a, ϵ_2_b
end    
    



################
# Floquet code #
################

# Define data form: ΔE_n | J | n_photons_a | n_photons_b
floquet_data_array = zeros( (N^2)*length(g_array)*length(g_4_array), 4 )

counter = 1
pbar = ProgressBar(total=length(g_4_array)*length(g_array))
for (_,g_4) in enumerate(g_4_array)
    N, T, θ, ω_0_a, ω_0_b, g_3_a, g_4_a, g_3_b, g_4_b, Ω_p, ω_p, Ω_d, ω_d, ϵ_2_a, ϵ_2_b = new_params(N,0,ω_0_a_og,ω_0_b_og,g_3_a_og,g_4,g_3_b_og,g_4) 

    # The Kerrs, needed to construct the dollowing coherent states
    K_a = (10*g_3_a^2)/(3*ω_0_a) - 3*g_4_a/2
    K_b = (10*g_3_b^2)/(3*ω_0_b) - 3*g_4_b/2

    # Propose a coherent superposition as initial groundstate of each well
    α_A = √(Complex(ϵ_2_a/K_a))
    ψ_0_a = coherent_state(N,α_A) + coherent_state(N,-α_A)
    ψ_0_a = ψ_0_a/norm(ψ_0_a)

    α_B = √(Complex(ϵ_2_b/K_b))
    ψ_0_b = coherent_state(N,α_B) + coherent_state(N,-α_B)
    ψ_0_b = ψ_0_b/norm(ψ_0_b)

    V_0 = kron(ψ_0_a,ψ_0_b)

    for (_,g) in enumerate(g_array)

        N, T, θ, ω_0_a, ω_0_b, g_3_a, g_4_a, g_3_b, g_4_b, Ω_p, ω_p, Ω_d, ω_d, ϵ_2_a, ϵ_2_b = new_params(N,g,ω_0_a_og,ω_0_b_og,g_3_a_og,g_4,g_3_b_og,g_4) 

        # display(new_params(N,g,ω_0_a_og,ω_0_b_og,g_3_a_og,g_4,g_3_b_og,g_4) )

        # find closest groundstate and eigenergies and update
        ϵ_n, ϕ_n = qen_qmodes_2_snails_dress(N,T, θ, ω_0_a, ω_0_b, g_3_a, g_4_a, g_3_b, g_4_b, Ω_p, ω_p, Ω_d, ω_d)

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
            ent = entanglement_entropy(ϕ_n[:,n])

            floquet_data_array[counter,:] .= Δϵ_n[n], g, g_4, ent
            counter += 1
        end
        update(pbar)
    end
end

# Put data in convenient DataFrame object and save it
df_floquet_2_snails = DataFrame(floquet_data_array, ["ΔE_n","g","g_4","ent"]) 
# CSV.write("data/floquet_crossings_comparison.csv", df_floquet)

# Combine both sims into a single dataframe
# df_comparison = DataFrame(vcat(floquet_data_array,H_eff_data_array), ["ΔE_n","ϵ_1","ϵ_2","n_photons","Floquet"]) 

# CSV.write("data/crossings_comparison.csv", df_comparison)
heatmap_matrix= zeros(length(g_array),length(g_4_array))
for (j,g) in enumerate(g_array)
    temp = filter(row ->  row.g == g, df_floquet_2_snails)
    for (k,g_4) in enumerate(g_4_array)
        temp_2 = filter(row ->  row.g_4 == g_4, temp)
        min_e = minimum(temp_2.ΔE_n)
        temp_3 = filter(row ->  row.ΔE_n == min_e, temp_2)
        heatmap_matrix[j,k] = temp_3.ent[1]
    end
end

plot!(log.(heatmap_matrix[:,1]))
for i=2:10
    plot!(log.( heatmap_matrix[:,i]))
    scatter!(log.( heatmap_matrix[:,i]))
end