include("../src/src.jl") # import src.jl which has creation/annahilation operators defined

using DataFrames # this is like pandas
using CSV 
using ProgressBars

# define parameter space in units of ω_0
N = 20
ω_0 = 1

g_n = [0.00075, 1.27*10^(-7)].*ω_0
Ω_1_array = Vector(0:0.00001:10^(-3)).*ω_0
Ω_2_array = Vector(0:0.05:6).*ω_0
ω_1 = ω_0*ω_0
ω_2 = 2*ω_0
cross_tol = ω_2/(2*1000) # distance at which two levels are considered to have crossed
n_cross = 3
# g_n = [2*10^(-5),8*10^(-6)]
# g_n = [1,1]


K = (10*g_n[1]^2)/(3*ω_0) - 3*g_n[2]/2

# # Define data form
# data_array = zeros( (N-1)*length(Ω_1_array)*length(Ω_2_array), 10 ) # data form: Δnn | ω_0 | Ω_1 | ω_1 | Ω_2 | ω_2 | K | ϵ_2 | ϵ_1 | N

# # The following is used to find ϵ_0 as described in the paper for your eyes only
# E_n , V_n = eigen(H_0(N,ω_0,g_n))
# V_0 = V_n[:, 1]

# counter = 1
# pbar = ProgressBar(total=length(Ω_1_array)*length(Ω_2_array))
# for Ω_1 in Ω_1_array
#     for Ω_2 in Ω_2_array

#         ϵ_n, ϕ_n = qen_qmodes(N, ω_0, g_n, Ω_1, ω_1, Ω_2, ω_2)

#         overlaps = zeros(N)
#         for n=1:N
#             overlaps[n] = real( (V_0'*ϕ_n[:,n])*(ϕ_n[:,n]'*V_0) )
#         end

#         ϵ_0 = ϵ_n[argmax(overlaps)]
#         V_0 = ϕ_n[:,argmax(overlaps)]
        
#         ϵ_1 = 2*Ω_1
#         ϵ_2 = 2*Ω_2*g_n[1]/(3*ω_0)
#         for n=1:(N-1)
#             Δnn = mod(  mod(ϵ_n[n+1] - ϵ_0, ω_2/2) - mod(ϵ_n[n] - ϵ_0, ω_2/2),  ω_2/2)
#             data_array[counter,:] .= Δnn, ω_0, Ω_1, ω_1,  Ω_2, ω_2, K, ϵ_2, ϵ_1, N
#             counter += 1
#         end
#         update(pbar)
#     end
# end


# Define data form
data_array = zeros( n_cross*length(Ω_1_array)*length(Ω_2_array), 10 ) # data form: Δnn | ω_0 | Ω_1 | ω_1 | Ω_2 | ω_2 | K | ϵ_2 | ϵ_1 | N

# The following is used to find ϵ_0 as described in the paper for your eyes only
E_n , V_n = eigen(H_0(N,ω_0,g_n))
V_0 = V_n[:, 1]

counter = 1
pbar = ProgressBar(total=length(Ω_1_array)*length(Ω_2_array))
for Ω_1 in Ω_1_array
    for Ω_2 in Ω_2_array

        ϵ_n, ϕ_n = qen_qmodes(N, ω_0, g_n, Ω_1, ω_1, Ω_2, ω_2)

        overlaps = zeros(N)
        for n=1:N
            overlaps[n] = real( (V_0'*ϕ_n[:,n])*(ϕ_n[:,n]'*V_0) )
        end

        ϵ_0 = ϵ_n[argmax(overlaps)]
        V_0 = ϕ_n[:,argmax(overlaps)]
        
        ϵ_1 = 2*Ω_1
        ϵ_2 = 2*Ω_2*g_n[1]/(3*ω_0)
        for n=1:n_cross
            Δnn = mod(  mod(ϵ_n[n+1] - ϵ_0, ω_2/2) - mod(ϵ_n[n] - ϵ_0, ω_2/2),  ω_2/2)
            data_array[counter,:] .= Δnn, ω_0, Ω_1, ω_1,  Ω_2, ω_2, K, ϵ_2, ϵ_1, N
            counter += 1
        end
        update(pbar)
    end
end



# Put data in convenient DataFrame object and save it
df_floquet = DataFrame(data_array, ["Δnn","ω_0","Ω_1","ω_1","Ω_2","ω_2","K","ϵ_2","ϵ_1","N"]) 
df_formatted = filter(row -> row.Δnn <= cross_tol, df_floquet)  # otherwise file is huge
CSV.write("data/floquet_resonances.csv", df_formatted)