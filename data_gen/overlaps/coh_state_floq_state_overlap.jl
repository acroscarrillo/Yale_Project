include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

using DataFrames # this is like pandas
using CSV 
using ProgressBars
using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs


N = 200

# in H_eff units:
Δ = 0
K = 1
ϵ_1_max = 70*K
ϵ_2_max = 30*K

# Define Floquet parameter space in units of ω_0
ω_0 = 1
ω_1 = 1
g_n = [0.00075, 1.27*10^(-7)].*ω_0
K = (10*g_n[1]^2)/(3*ω_0) - 3*g_n[2]/2
Ω_1_array = Vector( range(0, (1*K)*(ϵ_1_max)/2, length=150) )
Ω_2_array = Vector( range(0, 3*ϵ_2_max*(1*K)/(4*g_n[1]), length=200) )

# Define data form
data_array = zeros( n_cross*length(Ω_1_array)*length(Ω_2_array), 12 ) # data form: ϵ_n | Δnn | ω_0 | Ω_1 | ω_1 | Ω_2 | ω_2 | K | ϵ_2 | ϵ_1 | N | n_cross

counter = 1
pbar = ProgressBar(total=length(Ω_1_array)*length(Ω_2_array))
Ω_1, Ω_2 = 0, 0
for (i,_) in enumerate(Ω_1_array)
    for (j,_) in enumerate(Ω_2_array)
        if isodd(i)
            Ω_1, Ω_2 = Ω_1_array[i], Ω_2_array[j]
        else
            Ω_1, Ω_2 = Ω_1_array[i], Ω_2_array[length(Ω_2_array)-j+1]
        end

        ω_2 = 2*ω_a(ω_0,g_n,Ω_2)
    
        ϵ_2 = g_n[1]*Π(Ω_2,ω_2)
        ϵ_1 = 2*Ω_1

        T = ω_2/2
        E_n, v_n = H(N,ω_0,g_n,Ω_1,ω_1,Ω_2,ω_2,T)
        ipr_array = zeros(N)
        for n=1:N
            ipr_array[n] = IPR(v_n[:,n])
        end
    end
end