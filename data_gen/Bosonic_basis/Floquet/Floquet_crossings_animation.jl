include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

using LinearAlgebra # to diagonalise stuff
using DataFrames # this is like pandas
using ProgressBars
using CSV 
using Plots 


# Define parameter space
N = 20
Δ = 0
K = 1

ϵ_1_array = Vector(0.1:0.1:12)
ϵ_2_array = Vector(0:0.1:15)

n_crossings = 8

# Define crossing data form
data_array = zeros( n_crossings*length(ϵ_1_array)*length(ϵ_2_array), 7 ) # data form: ΔE_n | n |  Δ | K | ϵ_1 | ϵ_2 | N,  where Δnn is the difference between Nearest Neightbour levels


# Generate data within parameter space
counter = 1
for ϵ_1 in ProgressBar(ϵ_1_array)
    for ϵ_2 in ϵ_2_array
        H_temp = Hermitian( -H_eff(N,Δ,K,ϵ_1,ϵ_2) )
        lamb, _ = eigen(H_temp)
        for n=1:n_crossings
            data_array[counter,:] .= lamb[n]-minimum(lamb), n, Δ, K, ϵ_1, ϵ_2, N
            counter += 1
        end
    end
end

# put data in convenient DataFrame object
df_h_eff_repul = DataFrame(data_array, ["ΔE_n","n","Δ","K","ϵ_1","ϵ_2","N"]) 
CSV.write("data/h_eff_crossings_epsilon_1.csv", df_h_eff_repul)