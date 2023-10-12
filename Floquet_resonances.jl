#!/usr/bin/env julia

include("../src/src.jl") # import src.jl which has creation/annahilation operators defined

using DataFrames # this is like pandas
using CSV 

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

n_cross = N-1

# Define data form
data_array = zeros( n_cross*length(Ω_1_array)*length(Ω_2_array), 12 ) # data form: ϵ_n | Δnn | ω_0 | Ω_1 | ω_1 | Ω_2 | ω_2 | K | ϵ_2 | ϵ_1 | N | n_cross
df = DataFrame(data_array, ["ϵ_n","Δnn","ω_0","Ω_1","ω_1","Ω_2","ω_2","K","ϵ_2","ϵ_1","N","n_cross"]) 
CSV.write("testo.csv", df)

display("it worked!")