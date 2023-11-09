include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

using LinearAlgebra # to diagonalise stuff
using DataFrames # this is like pandas
using CSV 
using ProgressBars

N = 70 #it works just fine


##############
# H_eff code #
##############

# Define H_eff parameter space, in units of K
K = 11
Δ = -0*K
ϵ_1 = 0*K
ϵ_2_array = Vector(0:0.05:30).*K

# Define data form
H_eff_data_array = zeros( N*length(ϵ_2_array), 3 ) # data form: ΔE_n | ϵ_2 | Floquet?

# Generate H_eff data within parameter space
counter = 1
for ϵ_2 in ProgressBar(ϵ_2_array)
    H_temp = Hermitian( -H_eff(N,Δ,K,ϵ_1,ϵ_2) ) # H_eff/K as paper
    lamb, _ = eigen(H_temp)
    for i=1:N
        ΔE_n = lamb[i]-minimum(lamb)
        H_eff_data_array[counter,:] .= ΔE_n/K, ϵ_2/K, 0
        # H_eff_data_array[counter,:] .= ΔE_n, Δ, K, ϵ_1, ϵ_2, N, "H_eff"
        counter += 1
    end
end

# put data in convenient DataFrame object
# df_h_eff = DataFrame(H_eff_data_array, ["ΔE_n","Δ","K","ϵ_1","ϵ_2","N","theory"]) 
df_h_eff = DataFrame(H_eff_data_array, ["ΔE_n","ϵ_2","Floquet?"]) 
CSV.write("data/h_eff_kissing.csv", df_h_eff)


################
# Floquet code #
################

# Define Floquet parameter space
ω_0 = 1
# g_n = [0.00075, 1.27*10^(-7)].*ω_0 # exact match with H_eff
g_n = [-0.000975, -0.0000133].*ω_0 # experimental values
ω_1 = 0 .* ω_0
Ω_2_array = Vector(range(0,-0.7,length=200)).*ω_0 #graph above could be truncated at 10
Ω_1 = 0 .* ω_0
K = (10*g_n[1]^2)/(3*ω_0) - 3*g_n[2]/2


# Define data form
floquet_data_array = zeros( N*length(Ω_2_array), 3 ) # data form: Δϵ_n | ω_0 | g_3 | g_4 | Ω_1 | ω_1 | Ω_2 | ω_2 | K | ϵ_2/K | N | theory

# The following is used to find ϵ_0 as described in the paper for your eyes only
E_n , V_n = eigen(H_0(N,ω_0,g_n))
V_0 = V_n[:, 1]

# H_temp = Hermitian( H(N,ω_0,g_n,0,0,0) )
# lamb, _ = eigen(H_temp)
# ω_a = lamb[2] - lamb[1]
# ω_2 = 2*ω_a


counter = 1
for Ω_2 in ProgressBar(Ω_2_array) 
    ω_2 = 2*ω_a(ω_0,g_n,Ω_2)

    Π_temp = (Ω_2*ω_2)/(ω_2^2-ω_0^2)
    # display(abs(Π_temp/Π(Ω_2,ω_2))) # = 4.000something
    # Π_temp = Π(Ω_2,ω_2)
    ϵ_2 = g_n[1]*Π_temp

    ϵ_n, ϕ_n = qen_qmodes(N, ω_0, g_n, Ω_1, ω_1, Ω_2, ω_2)
    overlaps = zeros(N)
    for n=1:N
        overlaps[n] = real( (V_0'*ϕ_n[:,n])*(ϕ_n[:,n]'*V_0) )
    end

    ######################################################################
    # CAREFUL YOU DIVIDE BY 2K BUT IT SHOULD BE BY K. NEED TO CHECK THAT #
    ######################################################################

    ϵ_0 = ϵ_n[argmax(overlaps)]
    V_0 = ϕ_n[:,argmax(overlaps)]
    for n=1:N
        Δϵ_n = mod(ϵ_n[n] - ϵ_0, ω_2/2)
        floquet_data_array[counter,:] .= Δϵ_n/K, ϵ_2/K, 1

        # floquet_data_array[counter,:] .= Δϵ_n/K, ω_0, g_n[1], g_n[2], Ω_1, ω_1, Ω_2, ω_2, K, ϵ_2/K, N, "Floquet"
        counter += 1
    end
end

# Put data in convenient DataFrame object and save it
df_floquet = DataFrame(floquet_data_array, ["ΔE_n","ϵ_2","Floquet?"]) 
CSV.write("data/floquet_kissing.csv", df_floquet)


# Combine both sims into a single dataframe
df_comparison = DataFrame(vcat(floquet_data_array,H_eff_data_array), ["ΔE_n","ϵ_2","Floquet"]) 
CSV.write("data/kissing_comparison.csv", df_comparison)