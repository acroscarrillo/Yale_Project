include("../../src/src.jl") # import src.jl which has creation/annahilation operators, H_eff...etc defined

using DataFrames # this is like pandas
using CSV 
using ProgressBars


function H_quimico(N,k_1,k_2,k_4)
    A = a(N)
    x = A+A'
    p = im*(A-A')
    return p^2 + k_1*x -k_2*x^2 + k_4*x^4
end


# define parameter space
N = 200
k_4 = 1
k_1_array = Vector(0.1:0.05:14)
k_2_array = Vector(5:0.1:25)
anticross_n = 18

# Define crossing data form
crossing_data = zeros( anticross_n*length(k_1_array)*length(k_2_array), 7 ) # data form: cross_n |Δnn | Δ | k_4 | k_1 | k_2 | N,  where Δnn is the difference between Nearest Neightbour levels

# Generate data within parameter space
counter = 1 # Im lazy, sorry (it's just easier to read)
pbar = ProgressBar(total=length(k_1_array)*length(k_2_array))
for k_1 in k_1_array
    for k_2 in k_2_array
        H_temp = Hermitian( H_quimico(N,k_1,k_2,k_4) )
        lamb, _ = eigen(H_temp)
        for n=1:anticross_n # as anticrossings are: (1,2), (2,3), ...etc.
            crossing_data[counter,:] .= n, lamb[n+1] - lamb[n], Δ, k_4, k_1, k_2, N
            counter += 1
        end
        update(pbar)
    end
end

# put crossing data in convenient DataFrame object & save it
df_2_save = DataFrame(crossing_data, ["cross_n","Δnn", "Δ", "k_4", "k_1", "k_2","N"]) 

CSV.write("data/h_chem_anticrossings_lvl_tracking.csv", df_2_save)