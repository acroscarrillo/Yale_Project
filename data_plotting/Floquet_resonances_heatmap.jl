include("../src/src.jl") # import src.jl which has creation/annahilation operators defined

using DataFrames # this is like pandas
using CSV 
using ProgressBars
using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs

#generate the Ω_i arrays
N = 30

# in H_eff units:
Δ = 0
K = 1
ϵ_1_max = 70*K
ϵ_2_max = 32*K

# Define Floquet parameter space in units of ω_0
ω_0 = 1
g_n = [-0.0025, -6.667*10^(-5)].*ω_0
K = (10*g_n[1]^2)/(3*ω_0) - 3*g_n[2]/2
Ω_1_array = Vector( range(0, (1*K)*(ϵ_1_max)/2, length=150) )
Ω_2_array = Vector( range(0, 3*ϵ_2_max*(1*K)/(4*g_n[1]), length=200) )

Ω_1_max = Ω_1_array[end]
Ω_2_max = Ω_2_array[end]

ϵ_1_max = 2*Ω_1_max/K # check: in K units like in H_eff above (replace 2K -> K)
ϵ_2_max = Ω_2_max*4*g_n[1]/(3*ω_0*K) 

########

cross_cutoff = 5e-3
df_floquet = DataFrame(CSV.File("data/floquet_resonances_a_matrix_elements.csv"))
df = filter(row -> row.Δnn <= cross_cutoff, df_floquet)  # discard large crossings

n_cross_cutoff = 30
df = filter(row ->  row.n_cross <= n_cross_cutoff, df)  # discard large n's
df = filter(row ->  row.nn_cross <= n_cross_cutoff, df)  # discard large n's
# a_dag_nm_threshold = 1e-3
# df = filter(row ->  row.a_dag_nm < a_dag_nm_threshold, df_formatted)  # discard to-be-small couplings with bath

# photon_n_threshold = 5
# df = filter(row ->  row.photon_n < photon_n_threshold, df_formatted)   # discard states with large occupation number


matrix_color_2 = zeros(length(Ω_1_array),length(Ω_2_array))

pbar = ProgressBar(total=length(Ω_1_array)*length(Ω_2_array))
Ω_11, Ω_22 = 0, 0
for (i,_) in enumerate(Ω_1_array)
    for (j,_) in enumerate(Ω_2_array)
        if isodd(i)
            Ω_11, Ω_22 = Ω_1_array[i], Ω_2_array[j]
        else
            Ω_11, Ω_22 = Ω_1_array[i], Ω_2_array[length(Ω_2_array)-j+1]
        end

        min_temp = cross_cutoff
        try
            min_temp = minimum( filter(row ->  row.Ω_1 == Ω_11 && row.Ω_2 == Ω_22, df).Δnn )
        catch
            min_temp = cross_cutoff
        end
        Ω_1_ind = findfirst(Ω_1_array .== Ω_11)
        Ω_2_ind = findfirst(Ω_2_array .== Ω_22)
        matrix_color_2[Ω_1_ind,Ω_2_ind] = min_temp
        update(pbar)
    end
end

# ttl = "Floquet reson. at "*L"N="*string(N)*", "*L"\omega_0="*string(ω_0)*", "*L"\omega_1="*string(ω_1)*", "*L"\omega_2 = 2\omega_a"*", \n"*L"K="*string(round(K,sigdigits=3))*", "*L"g_n="*g_n*". All levels. Lin-scale."

ttl  = "N = $N , n_cross <= $n_cross_cutoff" #, a_dag_nm_threshold = $a_dag_nm_threshold, photon_n = $photon_n_threshold"

K = (10*g_n[1]^2)/(3*ω_0) - 3*g_n[2]/2
ϵ_1_array = Vector( range(0, ϵ_1_max*K, length=150) ) #convert back to ω_0 units
ϵ_2_array = Vector( range(0, ϵ_2_max*K, length=200) ) #convert back to ω_0 units
heatmap(ϵ_1_array,ϵ_2_array, log.(matrix_color_2'), xlab=L"\Omega_1",ylab=L"\Omega_2",guidefontsize=14, title = ttl)#,size=(700,600))


df_floquet = DataFrame(CSV.File("data/floquet_resonances_full_backup.csv"))
df_Omega = filter(row -> row.Ω_2 ==  0.014474849246231157, df_floquet)  # discard large crossings

df_Omega_2 = filter(row -> row.n_cross <70 && row.Δnn <= 0.00005, df_Omega)  # discard large crossings

df_Omega_2 = filter(row -> row.n_cross <70 && row.ϵ_n <= 0.001, df_Omega)  # discard large crossings

df_Omega_2 = filter(row -> row.n_cross <5, df_Omega)  # discard large crossings


@df df_Omega_2 plot(:ϵ_1, :ϵ_n, group=:n_cross, markersize = 1, xlab = L"\epsilon_1",ylab = L"\epsilon_n",seriestype = :scatter,alpha=1,guidefontsize=14,markerstrokewidth=0, legend=true,legendtitle = "Levels",foreground_color_legend = nothing,background_color_legend=nothing)