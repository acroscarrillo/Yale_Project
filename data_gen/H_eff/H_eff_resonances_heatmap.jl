include("../../src/src.jl") # import src.jl which has creation/annahilation operators, H_eff...etc defined

using DataFrames # this is like pandas
using CSV 
using ProgressBars

# define parameter space
N = 200
Δ = 0
K = 1
ϵ_1_array = Vector(0.1:0.03:12)
ϵ_2_array = Vector(0:0.07:10)
cross_tol = 0.2 # distance at which two levels are considered to have crossed

# Define crossing data form
crossing_data = zeros( length(ϵ_1_array)*length(ϵ_2_array), 6 ) # data form: min(Δnn) | Δ | K | ϵ_1 | ϵ_2 | N ,  where Δnn is the difference between Nearest Neightbour levels

matrix_color = zeros(length(ϵ_2_array),length(ϵ_1_array))

# Generate data within parameter space
counter = 1 # Im lazy, sorry (it's just easier to read)
pbar = ProgressBar(total=length(ϵ_1_array)*length(ϵ_2_array))
for (i,ϵ_1) in enumerate(ϵ_1_array)
    for (j,ϵ_2) in enumerate(ϵ_2_array)
        H_temp = Hermitian( -H_eff(N,Δ,K,ϵ_1,ϵ_2) )
        lamb, _ = eigen(H_temp)
        temp_cross =  zeros(N-1)
        for n=1:(N-1)
            temp_cross[n] = lamb[1+n] - lamb[n]
        end
        crossing_data[counter,:] .= minimum(temp_cross), Δ, K, ϵ_1, ϵ_2, N
        matrix_color[j,i] = minimum(temp_cross)
        counter += 1
        update(pbar)
    end
end

# put crossing data in convenient DataFrame object & save it
df_crossing = DataFrame(crossing_data, ["min(Δnn)", "Δ", "K", "ϵ_1", "ϵ_2","N"]) 
CSV.write("data/h_eff_heatmap_crossing.csv", df_crossing)