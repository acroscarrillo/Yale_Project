include("../../src/src.jl") # import src.jl which has creation/annahilation operators, H_eff...etc defined

using DataFrames # this is like pandas
using CSV 
using ProgressBars


# define parameter space
N = 200
Δ = 0
K = 1
ϵ_1_array = Vector(range(0.1, 13.7, length=100))
ϵ_2_array = Vector(range(2.5, 13, length=100))

# define FG rule conditions
perturb = a(N)'
i = 3
f = 1

# simulation
mat_elems_array = zeros(length(ϵ_1_array),length(ϵ_2_array))
pbar = ProgressBar(total=length(ϵ_1_array)*length(ϵ_2_array))
for (j,ϵ_1) in enumerate(ϵ_1_array)
    for (k,ϵ_2) in enumerate(ϵ_2_array)
        H_temp = H_eff(N,Δ,K,ϵ_1,ϵ_2)
        E, V = eigen(H_temp)
        mat_elems_array[j,k] = norm(V[:,i]'*perturb*V[:,f])
        update(pbar)
    end
end

heatmap(mat_elems_array)