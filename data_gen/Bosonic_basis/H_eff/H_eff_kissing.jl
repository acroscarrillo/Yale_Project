include("../src/src.jl") # import src.jl which has creation/annahilation operators defined

using LinearAlgebra # to diagonalise stuff
using DataFrames # this is like pandas
using CSV 

function H_qu(N, ξ)
    A = a(N) # annahilation op. up to dim N
    return  (A'^2)*(A^2) -  ξ*(A'^2 + A^2) 
end


# Define parameter space
N = 100
ξ_array = Vector(0:0.01:12)

# Define data form
data_array = zeros( N*length(ξ_array), 3 ) # data form: ΔE_n | ξ | N

# Generate data within parameter space
counter = 1
for ξ in ξ_array
    H_temp = Hermitian( H_qu(N, ξ) )
    lamb, _ = eigen(H_temp)
    for i=1:N
        data_array[counter,:] .= lamb[i]-minimum(lamb), ξ, N
        counter += 1
    end
end

# put data in convenient DataFrame object
df_h_eff = DataFrame(data_array, ["ΔE_n","ξ","N"]) 
CSV.write("data/h_eff_kissing.csv", df_h_eff)

