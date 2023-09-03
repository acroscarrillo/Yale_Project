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

# Define H_sk model as in https://arxiv.org/pdf/2209.03934.pdf
function H_sk(N, ϵ_2_by_K)
    return H_eff(N, 0, 1, 0, ϵ_2_by_K)
end

# Define parameter space
N = 20
ϵ_2_by_K_array = Vector(0:0.2:40)

# Define data form
data_array = zeros( N*length(ϵ_2_by_K_array), 3 ) # E_n  ϵ_2_by_K  N

# Generate data within parameter space
counter = 1
for ϵ_2_by_K in ProgressBar(ϵ_2_by_K_array)
    H = Hermitian( H_sk(N, ϵ_2_by_K) )
    lamb, vec_array = eigen(H)
    for i=1:N
        data_array[counter,:] .= lamb[i], ϵ_2_by_K, N
        counter += 1
    end
end

# Store data
df_wb = DataFrame(data_array, ["E_n","ϵ_2_by_K","N"])
CSV.write("data/H_sk_spectrum.csv", df_wb)
