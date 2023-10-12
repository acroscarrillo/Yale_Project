include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

using DataFrames # this is like pandas
using CSV 
using ProgressBars

N = 30    

# in H_eff units: 
Δ = 0
K = 1
ϵ_1_max = 70*K
ϵ_2_max = 30*K

# Define Floquet parameter space in units of ω_0

ω_0 = 1
# ω_1 = 1
# g_n = [0.00075, 1.27*10^(-7)].*ω_0
g_n = [-0.0025, -6.667*10^(-5)].*ω_0
K = (10*g_n[1]^2)/(3*ω_0) - 3*g_n[2]/2
Ω_1_array = Vector( range(0, (1*K)*(ϵ_1_max)/2, length=150) )
Ω_2_array = Vector( range(0, 3*ϵ_2_max*(1*K)/(4*g_n[1]), length=200) )

Ω_1_max = Ω_1_array[end]
Ω_2_max = Ω_2_array[end]

ϵ_1_max = 2*Ω_1_max/K # check: in K units like in H_eff above (replace 2K -> K)
ϵ_2_max = Ω_2_max*4*g_n[1]/(3*ω_0*K) # check: in K units like in H_eff above (replace 2K -> K)

# cross_tol = ω_2/(2*100) # distance at which two levels are considered to have crossed
n_cross = N-1


# Define data form
data_array = zeros( n_cross*length(Ω_1_array)*length(Ω_2_array), 16 ) # data form: ϵ_n | Δnn | ω_0 | Ω_1 | ω_1 | Ω_2 | ω_2 | K | ϵ_2 | ϵ_1 | N | n_cross | a_nm

# The following is used to find ϵ_0 as described in the paper for your eyes only
E_0 , V_0_array = eigen(H_0(N,ω_0,g_n))


# H_temp = Hermitian( H(N,ω_0,g_n,0,0,0) )
# lamb, _ = eigen(H_temp)
# ω_a = lamb[2] - lamb[1]
# ω_2 = 2*ω_a

matrix_color = zeros(length(Ω_1_array),length(Ω_2_array))

counter = 1
pbar = ProgressBar(total=length(Ω_1_array)*length(Ω_2_array))
Ω_1, Ω_2 = 0, 0
for (i,_) in enumerate(Ω_1_array)
    for (j,_) in enumerate(Ω_2_array)
        if isodd(i)
            Ω_1, Ω_2 = Ω_1_array[i], Ω_2_array[j]
        else
            Ω_1, Ω_2 = Ω_1_array[i], Ω_2_array[length(Ω_2_array)-j+1]
        end

        ω_2 = 2*ω_a(ω_0,g_n,Ω_2)
        ω_1 = ω_2/2
        ϵ_2 = g_n[1]*Π(Ω_2,ω_2)
        ϵ_1 = 2*Ω_1
        
        ϵ_n, ϕ_n = qen_qmodes(N, ω_0, g_n, Ω_1, ω_1, Ω_2, ω_2)
        overlaps_array = norm.( (V_0_array'*ϕ_n) ) # get matrix of overlaps
        V_0_array = ϕ_n # update old vectors

        levels_order_array = mapslices(argmax,overlaps_array,dims=2) # vector of levels positions, i.e. the groundstate position is levels_order_array[1].

        ϵ_0 = ϵ_n[levels_order_array[1]]

        ϵ_n = sort( mod.( ϵ_n .- ϵ_0, ω_2/2) ) # sort them all properly
        for n=1:n_cross
            Δnn = mod( ϵ_n[n+1] - ϵ_n[n],  ω_2/2)
            a_dag_nm = norm( ϕ_n[:,n+1]'*a(N)'*ϕ_n[:,n] )
            a_nm = norm( ϕ_n[:,n+1]'*a(N)*ϕ_n[:,n] )
            photon_n = real( ϕ_n[:,n]'*a(N)'*a(N)*ϕ_n[:,n] )
            data_array[counter,:] .= ϵ_n[n+1], Δnn, ω_0, Ω_1, ω_1, Ω_2, ω_2, K, ϵ_2, ϵ_1, N, levels_order_array[n+1],levels_order_array[n], a_dag_nm, a_nm, photon_n
            counter += 1
        end

        temp_cross = zeros(n_cross)
        for n=1:(n_cross)
            temp_cross[n] =  mod( ϵ_n[n+1] - ϵ_n[n],  ω_2/2)
        end
        Ω_1_ind = findfirst(Ω_1_array .== Ω_1)
        Ω_2_ind = findfirst(Ω_2_array .== Ω_2)
        matrix_color[Ω_1_ind,Ω_2_ind] = minimum(temp_cross)

        update(pbar)
    end
end

K = (10*g_n[1]^2)/(3*ω_0) - 3*g_n[2]/2
ϵ_1_array = Vector( range(0, ϵ_1_max*K, length=150) ) #convert back to ω_0 units
ϵ_2_array = Vector( range(0, ϵ_2_max*K, length=200) ) #convert back to ω_0 units
heatmap(ϵ_1_array,ϵ_2_array,log.(matrix_color'), xlab=L"\epsilon_1/\omega_0",ylab=L"\epsilon_2/\omega_0",title=L"N="*string(N))

# ϵ_1_array = Vector( range(0, ϵ_1_max, length=150) )
# ϵ_2_array = Vector( range(0, ϵ_2_max, length=200) )
# heatmap(ϵ_1_array,ϵ_2_array,log.(matrix_color'), xlab=L"\epsilon_1/K",ylab=L"\epsilon_2/K",title=L"N="*string(N))

# Put data in convenient DataFrame object and save it
df_floquet = DataFrame(data_array, ["ϵ_n","Δnn","ω_0","Ω_1","ω_1","Ω_2","ω_2","K","ϵ_2","ϵ_1","N","nn_cross","n_cross","a_dag_nm","a_nm","photon_n"]) 
# df_formatted = filter(row -> row.Δnn <= 0.0000032, df_floquet)  # discard large crossings
# df_formatted = filter(row -> row.Δnn <= 0.000005, df_floquet)  # discard large crossings

CSV.write("data/floquet_resonances.csv", df_formatted)

# df_formatted = filter(row ->  row.n_cross == 1 && row.Ω_1 == 0 && row.Ω_2 == 0, df_floquet)