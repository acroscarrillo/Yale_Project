include("../../src/main_charge_basis.jl") # import charge basis stuff
using ProgressBars
using DataFrames

N = 101
E_J = 1  
E_C = 1/5907
m = 2
α = 0.1
ϕ_ext = 0.33*2*π


###########
# Floquet #
###########
# Define parameter space
ω_1 = 1
ω_2 = 2
Ω_1_array = Vector( range(0, 0, length=2) )
Ω_2_array = Vector( range(0, 1e-1, length=200) )

# define data array to store Floquet data
floquet_array = zeros(N*length(Ω_1_array)*length(Ω_2_array), 4 )  # data form: ΔE_n | ϵ_2 | ϵ_1 | floquet?

# Generate data within parameter space
E_n , V_n  = eigen(H_joseph(N,E_J,α,ϕ_ext,m))
V_0 = V_n[:, 1]

heatmap_mat = zeros(length(Ω_2_array),length(Ω_1_array))

pbar = ProgressBar(total=length(Ω_1_array)*length(Ω_2_array))
counter = 1
Ω_1, Ω_2 = 0, 0
for (i,_) in enumerate(Ω_2_array)
    for (j,_) in enumerate(Ω_1_array)
        # otherwise for loop is not adiabatic (makes jumps in Ω_1, Ω_2), this is like a snake in the Ω_1-Ω_2 plane. Failing to do this is an error.
        if isodd(i)
            Ω_2, Ω_1 = Ω_2_array[i], Ω_1_array[j]
        else
            Ω_2, Ω_1 = Ω_2_array[i], Ω_1_array[length(Ω_1_array)-j+1]
        end

        ϵ_n, ϕ_n = qen_qmodes(N, E_C, E_J, α, ϕ_ext, m, Ω_1, ω_1, Ω_2, ω_2)

        # adiabatic tracking of the ground state
        overlaps = zeros(N)
        for n=1:N
            overlaps[n] = real( (V_0'*ϕ_n[:,n])*(ϕ_n[:,n]'*V_0) )
        end

        ϵ_0 = ϵ_n[argmax(overlaps)]
        V_0 = ϕ_n[:,argmax(overlaps)]
        ϵ_n = sort(mod.(ϵ_n .- ϵ_0, ω_2/2))

        for n=1:N
            floquet_array[counter,:] .= ϵ_n[n]/(1*E_J), Ω_1/(1*E_J), Ω_2/(1*E_J), 1
            counter += 1
        end
        heatmap_mat[i,j] = minimum(ϵ_n)
        update(pbar)
    end
end


# put data in convenient DataFrame object
df_floquet = DataFrame(floquet_array, ["ΔE_n","Ω_1","Ω_2","floquet"]) 
# CSV.write("data/charge_basis.csv", df_floquet)

heatmap(heatmap_mat)

df_temp = filter(row ->  row.Ω_1 == 0, df_floquet)
scatter(df_temp.Ω_2,df_temp.ΔE_n,ms=0.5,ylim=(0,1))