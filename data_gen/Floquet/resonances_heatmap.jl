include("../../src/src.jl") # import src.jl which has creation/annahilation operators defined

using DataFrames # this is like pandas
using CSV 
using ProgressBars
using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs

#generate the Ω_i arrays
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
########


df_floquet = DataFrame(CSV.File("data/floquet_resonances_full_backup.csv"))
df_formatted = filter(row -> row.Δnn <= 0.005, df_floquet)  # discard large crossings

heat_cross = 200
df = filter(row ->  row.n_cross < heat_cross, df_formatted)  # match paper limits

matrix_color_2 = ones(length(Ω_1_array),length(Ω_2_array))

pbar = ProgressBar(total=length(Ω_1_array)*length(Ω_2_array))
Ω_11, Ω_22 = 0, 0
for (i,_) in enumerate(Ω_1_array)
    for (j,_) in enumerate(Ω_2_array)
        if isodd(i)
            Ω_11, Ω_22 = Ω_1_array[i], Ω_2_array[j]
        else
            Ω_11, Ω_22 = Ω_1_array[i], Ω_2_array[length(Ω_2_array)-j+1]
        end

        min_temp = 0.005
        try
            min_temp = minimum( filter(row ->  row.Ω_1 == Ω_11 && row.Ω_2 == Ω_22, df).Δnn )
        catch
            min_temp = 0.005
        end
        Ω_1_ind = findfirst(Ω_1_array .== Ω_11)
        Ω_2_ind = findfirst(Ω_2_array .== Ω_22)
        matrix_color_2[Ω_1_ind,Ω_2_ind] = min_temp
        update(pbar)
    end
end

heatmap(Ω_1_array,Ω_2_array, matrix_color_2', xlab=L"\Omega_1",ylab=L"\Omega_2",title="N=200, n_levels="*string(heat_cross))


df_floquet = DataFrame(CSV.File("data/floquet_resonances_full_backup.csv"))
df_Omega = filter(row -> row.Ω_2 ==  0.035298316582914574, df_floquet)  # discard large crossings

df_Omega_2 = filter(row -> row.n_cross <70 && row.Δnn <= 0.00005, df_Omega)  # discard large crossings

df_Omega_2 = filter(row -> row.n_cross <70 && row.ϵ_n <= 0.014  && row.ϵ_n >= 0.012, df_Omega)  # discard large crossings

df_Omega_2 = filter(row -> row.n_cross <70 && row.ϵ_n <= 0.014  && row.ϵ_n >= 0.012, df_Omega)  # discard large crossings


@df df_Omega_2 plot(:ϵ_1, :ϵ_n, group=:n_cross, markersize = 1, xlab = L"\epsilon_1",ylab = L"\epsilon_n",seriestype = :scatter,alpha=1,guidefontsize=14,markerstrokewidth=0, legend=true,legendtitle = "Levels",foreground_color_legend = nothing,background_color_legend=nothing)