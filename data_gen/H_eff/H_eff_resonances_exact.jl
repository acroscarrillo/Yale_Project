include("../../src/src.jl") # import src.jl which has creation/annahilation operators, H_eff...etc defined

using DataFrames # this is like pandas
using CSV 
using ProgressBars

# define parameter space
N = 200
Δ = 0
K = 1
ϵ_1_array = Vector(0:0.03:10)
ϵ_2_array = Vector(0:0.07:10)
cross_tol = 0.5 # distance at which two levels are considered to have crossed
n_crossings = 199

# # Define crossing data form
# crossing_data = zeros( length(ϵ_1_array)*length(ϵ_2_array), 6 ) # data form: min(Δnn) | Δ | K | ϵ_1 | ϵ_2 | N,  where Δnn is the difference between Nearest Neightbour levels


# # Generate data within parameter space
# counter = 1 # Im lazy, sorry (it's just easier to read)
# for ϵ_1 in ProgressBar(ϵ_1_array)
#     for ϵ_2 in ϵ_2_array
#         H_temp = Hermitian( H_eff(N,Δ,K,ϵ_1,ϵ_2) )
#         lamb, _ = eigen(H_temp)
#         temp_cross =  zeros(N-1)
#         for n=1:(N-1)
#             temp_cross[n] = lamb[1+n] - lamb[n]
#         end
#         crossing_data[counter,:] .= minimum(temp_cross), Δ, K, ϵ_1, ϵ_2, N
#         counter += 1
#     end
# end


# Define crossing data form
crossing_data = zeros( n_crossings*length(ϵ_1_array)*length(ϵ_2_array), 7 ) # data form: cross_n |Δnn | Δ | K | ϵ_1 | ϵ_2 | N,  where Δnn is the difference between Nearest Neightbour levels


# Generate data within parameter space
counter = 1 # Im lazy, sorry (it's just easier to read)
pbar = ProgressBar(total=length(ϵ_1_array)*length(ϵ_2_array))
for ϵ_1 in ϵ_1_array
    for ϵ_2 in ϵ_2_array
        H_temp = Hermitian( -H_eff(N,Δ,K,ϵ_1,ϵ_2) )
        lamb, _ = eigen(H_temp)
        for n=1:n_crossings
            crossing_data[counter,:] .= n, lamb[1+n] - lamb[n], Δ, K, ϵ_1, ϵ_2, N
            counter += 1
        end
        update(pbar)
    end
end

temp_theory = zeros(+length(ϵ_1_array),7)
counter = 1
for ϵ_1 in ϵ_1_array
    ϵ_2 = (4*K*ϵ_1^2)^(1/3)
    temp_theory[counter,:] .= -1, ϵ_2, Δ, K, ϵ_1, ϵ_2, N
    counter += 1
end


# put crossing data in convenient DataFrame object & save it
df_crossing = DataFrame(crossing_data, ["cross_n","Δnn", "Δ", "K", "ϵ_1", "ϵ_2","N"]) 
df_formatted = filter(row -> row.Δnn <= cross_tol, df_crossing)  # otherwise file is huge

df_theory = DataFrame(temp_theory, ["cross_n","Δnn", "Δ", "K", "ϵ_1", "ϵ_2","N"]) 
df_formatted = vcat(df_formatted, df_theory)

CSV.write("data/h_eff_exact_resonances.csv", df_formatted)