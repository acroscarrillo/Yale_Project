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

# Define H_sk model as in (2) of https://arxiv.org/pdf/2209.03934.pdf
function H_sk(N, ϵ_2_by_K)
    A = a(N)
    return   ϵ_2_by_K*(A'^2 + A^2) + (A'^2)*(A^2) + (A' + A)
end

function H_sk_2(N, ϵ_2_by_K)
    A = a(N)
    n = a(N)'*a(N)
    return 2*n+I(N)-0.25*(2*n+I(N))^2 + ϵ_2_by_K*(A^2 + A'^2)
end

# Define parameter space
N = 40
ϵ_2_by_K_array = Vector(0:0.1:40)

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
CSV.write("data/H_qu_spectrum.csv", df_wb)
