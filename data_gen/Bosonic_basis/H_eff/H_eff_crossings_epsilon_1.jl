include("../src/src.jl") # import src.jl which has creation/annahilation operators defined

using LinearAlgebra # to diagonalise stuff
using DataFrames # this is like pandas
using ProgressBars
using CSV 

# Define parameter space
N = 30
Δ = 0
K = 1
ϵ_1_array = Vector(0:0.01:18)
ϵ_2 = 1


# Define crossing data form
data_array = zeros( N*length(ϵ_1_array), 6 ) # data form: ΔE_n | Δ | K | ϵ_1 | ϵ_2 | N,  where Δnn is the difference between Nearest Neightbour levels


# Generate data within parameter space
counter = 1
for ϵ_1 in ProgressBar(ϵ_1_array)
    H_temp = Hermitian( -H_eff(N,Δ,K,ϵ_1,ϵ_2) )
    lamb, _ = eigen(H_temp)
    for i=1:N
        data_array[counter,:] .= lamb[i]-minimum(lamb), Δ, K, ϵ_1, ϵ_2, N
        counter += 1
    end
end

# put data in convenient DataFrame object
df_h_eff_repul = DataFrame(data_array, ["ΔE_n","Δ","K","ϵ_1","ϵ_2","N"]) 
CSV.write("data/h_eff_crossings_epsilon_1.csv", df_h_eff_repul)

