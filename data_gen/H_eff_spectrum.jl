include("../src/src.jl")

using ProgressBars
using Statistics
using CSV 
using LinearAlgebra
using DataFrames

# Define H_eff model
function H_eff(N, Δ, K, ϵ_1, ϵ_2)
    A = a(N)
    return Δ*(A'*A) - K*(A'^2)*(A^2) + ϵ_1*(A' + A) + ϵ_2*(A'^2 + A^2)
end

# Define parameter space
N = 10
Δ_array = Vector(0:0)
K_array = Vector(-2.5:0.2:2.5)
ϵ_1_array = Vector(0:0.2:5)
ϵ_2_array = Vector(0:0.2:10)

# Define data form
data_array = zeros( N*length(Δ_array)*length(K_array)*length(ϵ_1_array)*length(ϵ_2_array), 6 ) # E_n  Δ  K  ϵ_1  ϵ_2  N

# Generate data within parameter space
counter = 1
for K in ProgressBar(K_array)
    for Δ in Δ_array
        for ϵ_1 in ϵ_1_array
            for ϵ_2 in ϵ_2_array
                H = Hermitian( H_eff(N, Δ, K, ϵ_1, ϵ_2) )
                lamb, vec_array = eigen(H)
                for i=1:N
                    data_array[counter,:] .= lamb[i], Δ, K, ϵ_1, ϵ_2, N
                    counter += 1
                end
            end
        end
    end
end

# Store data
df_wb = DataFrame(data_array, ["E_n","Δ","K","ϵ_1","ϵ_2","N"])
CSV.write("data/H_eff_spectrum.csv", df_wb)
