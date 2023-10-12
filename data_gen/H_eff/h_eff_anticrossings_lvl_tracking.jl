include("../../src/src.jl") # import src.jl which has creation/annahilation operators, H_eff...etc defined

using DataFrames # this is like pandas
using CSV 
using ProgressBars

# define parameter space
N = 200
Δ = 0
K = 1

ϵ_1_array = Vector(range(0.1, 13.7, length=300))
ϵ_2_array = Vector(range(2.5, 13, length=300))

anticross_n = 20

# Define crossing data form
crossing_data = zeros( anticross_n*length(ϵ_1_array)*length(ϵ_2_array), 7 ) # data form: cross_n |Δnn | Δ | K | ϵ_1 | ϵ_2 | N,  where Δnn is the difference between Nearest Neightbour levels

# Generate data within parameter space
counter = 1 # Im lazy, sorry (it's just easier to read)
pbar = ProgressBar(total=length(ϵ_1_array)*length(ϵ_2_array))
for ϵ_1 in ϵ_1_array
    for ϵ_2 in ϵ_2_array
        H_temp = Hermitian( -H_eff(N,Δ,K,ϵ_1,ϵ_2) )
        lamb, _ = eigen(H_temp)
        for n=1:anticross_n # as anticrossings are: (1,2), (2,3), ...etc.
            crossing_data[counter,:] .= n, lamb[n+1] - lamb[n], Δ, K, ϵ_1, ϵ_2, N
            counter += 1
        end
        update(pbar)
    end
end

# put crossing data in convenient DataFrame object & save it
df_2_save = DataFrame(crossing_data, ["cross_n","Δnn", "Δ", "K", "ϵ_1", "ϵ_2","N"]) 

CSV.write("data/data_to_gitignore/h_eff_anticrossings_lvl_tracking.csv", df_2_save)