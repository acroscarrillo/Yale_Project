using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs
using Plots

#########
# H_eff #
#########

df_H_eff = DataFrame(CSV.File("data/h_eff_heatmap_crossing.csv"))


N_eff = unique(df_H_eff.N)[1]
Δ_eff = unique(df_H_eff.Δ)[1]
K_eff = unique(df_H_eff.K)[1]

ϵ_1_array_eff = unique(df_H_eff.ϵ_1)
ϵ_2_array_eff = unique(df_H_eff.ϵ_2)

matrix_color_eff = reshape(df."min(Δnn)",(length(ϵ_2_array),length(ϵ_1_array)))

ttl_eff = L"H_{eff}"*" resonances at "*L"N="*string(N_eff)*", "*L"\Delta="*string(Δ_eff)*", "*L"K="*string(K_eff)*".\n Color bar in log scale. All levels considered."

h_eff_heatmap_plot = heatmap(ϵ_1_array_eff,ϵ_2_array_eff, log.(matrix_color_eff), xlab=L"\epsilon_1",ylab=L"\epsilon_2",title = ttl_eff, guidefontsize=9,titlefontsize=7)






h_eff_crossings_plot = plot(rand(100,3),layout = (1,3), axis=nothing,legend=false)
h_eff_waves_plot = plot(rand(100,3),layout = (1,3) , axis=nothing,legend=false)
plot(h_eff_heatmap_plot, h_eff_heatmap_plot, h_eff_waves_plot, 
h_eff_waves_plot, h_eff_crossings_plot, h_eff_waves_plot, 
layout = grid(3, 2, heights=[0.4 ,0.1, 0.4, 0.1]))