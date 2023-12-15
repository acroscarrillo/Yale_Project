include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

using LinearAlgebra # to diagonalise stuff
using DataFrames # this is like pandas
using CSV 
using ProgressBars


N = 150 #it works just fine for tests


##############
# H_eff code #
##############

# Define H_eff parameter space, in units of K
K = 1
Δ = -0*K
ϵ_1_array = Vector(range(0,14,length=200)).*K
ϵ_2_array = Vector(range(0,18,length=200)).*K
ϵ_1_array = ϵ_1_array[2:end]
ϵ_2_array = ϵ_2_array[2:end]

# Define data form: ΔE_n | ϵ_1 | ϵ_2 | Floquet?
H_eff_data_array = zeros( N*length(ϵ_1_array)*length(ϵ_2_array), 4 ) 

# Generate H_eff data within parameter space
counter = 1
for ϵ_1 in ProgressBar(ϵ_1_array)
    for ϵ_2 in ϵ_2_array
        H_temp = Hermitian( -H_eff(N,Δ,K,ϵ_1,ϵ_2) )
        lamb, _ = eigen(H_temp)
        for i=1:N
            ΔE_n = lamb[i]-minimum(lamb)
            H_eff_data_array[counter,:] .= ΔE_n/K, ϵ_1/K, ϵ_2/K, 0
            counter += 1
        end
    end
end

# put data in convenient DataFrame object
df_h_eff = DataFrame(H_eff_data_array, ["ΔE_n","ϵ_1","ϵ_2","Floquet?"]) 
CSV.write("data/h_eff_crossings_comparison.csv", df_h_eff)


################
# Floquet code #
################

# Define Floquet parameter space
ω_0 = 1
# g_n = [0.00075, 1.27*10^(-7)].*ω_0 # exact match with H_eff (garcia mata fig 1)
# g_n = [-0.0035, -3.33e-5].*ω_0 #
g_n = [-0.00331793, -5.0892e-5].*ω_0 # experimental values
K = (10*g_n[1]^2)/(3*ω_0) - 3*g_n[2]/2

Ω_1_array = Vector(range(0,Ω_1_max_by_ω_0(K,ϵ_1_array),length=200)).*ω_0 
Ω_2_array = Vector(range(0,Ω_2_max_by_ω_0(K,g_n,ϵ_2_array),length=200)).*ω_0
Ω_1_array = Ω_1_array[2:end]
Ω_2_array = Ω_2_array[2:end]


# Define data form: ΔE_n | ϵ_1 | ϵ_2 | Floquet?
floquet_data_array = zeros( N*length(Ω_1_array)*length(Ω_2_array), 4 )

# The following is used to find ϵ_0 adiabatically as in the paper for your eyes only
E_n , V_n = eigen(H_0(N,ω_0,g_n))
V_0 = V_n[:, 1]

counter = 1
pbar = ProgressBar(total=length(Ω_1_array)*length(Ω_2_array))
Ω_1, Ω_2 = 0, 0
for (i,_) in enumerate(Ω_1_array)
    for (j,_) in enumerate(Ω_2_array)
        # define adiabatic parameter path
        if isodd(i)
            Ω_1, Ω_2 = Ω_1_array[i], Ω_2_array[j]
        else
            Ω_1, Ω_2 = Ω_1_array[i], Ω_2_array[length(Ω_2_array)-j+1]
        end

        # define parameters
        ω_1, ω_2 = ω_a(ω_0,g_n,Ω_2), 2*ω_a(ω_0,g_n,Ω_2)
        ϵ_1, ϵ_2 = Ω_1/2, g_n[1]*Π(ω_0,Ω_2,ω_2)

        # find closest groundstate and eigenergies and update
        ϵ_n, ϕ_n = qen_qmodes(N, ω_0, g_n, Ω_1, ω_1, Ω_2, ω_2)
        # overlaps = zeros(N)
        # for n=1:N
        #     overlaps[n] = real( (V_0'*ϕ_n[:,n])*(ϕ_n[:,n]'*V_0) )
        # end
        # ϵ_0 = ϵ_n[argmax(overlaps)]
        # V_0 = ϕ_n[:,argmax(overlaps)]

        x_well, p_well = dblwell_minimas(K,ϵ_1,ϵ_2)[1], 0
        coh_stt = coherent_state(N,x_well,p_well)

        overlaps_array = norm.( (coh_stt'*ϕ_n) ) # get matrix of overlaps

        levels_order_array = mapslices(argmax,overlaps_array,dims=2) # vector of levels positions, i.e. the groundstate position is levels_order_array[1].

        ϕ_0 = ϕ_n[:,levels_order_array[1]]
        ϵ_0 = ϵ_n[levels_order_array[1]]

        Δϵ_n = sort( mod.( ϵ_n .- ϵ_0, ω_2/2) )
        #store parameters in array
        for n=1:N
            floquet_data_array[counter,:] .= Δϵ_n[n]/K, ϵ_1/K, ϵ_2/K, 1
            counter += 1
        end
        update(pbar)
    end
end

# Put data in convenient DataFrame object and save it
df_floquet = DataFrame(floquet_data_array, ["ΔE_n","ϵ_1","ϵ_2","Floquet?"]) 
CSV.write("data/floquet_crossings_comparison.csv", df_floquet)

# Combine both sims into a single dataframe
df_comparison = DataFrame(vcat(floquet_data_array,H_eff_data_array), ["ΔE_n","ϵ_1","ϵ_2","Floquet"]) 

CSV.write("data/crossings_comparison.csv", df_comparison)