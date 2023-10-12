using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs

df = DataFrame(CSV.File("data/h_eff_heatmap_crossing.csv"))

N = unique(df.N)[1]
Δ = unique(df.Δ)[1]
K = unique(df.K)[1]

ϵ_1_array = unique(df.ϵ_1)
ϵ_2_array = unique(df.ϵ_2)

matrix_color = reshape(df."min(Δnn)",(length(ϵ_2_array),length(ϵ_1_array)))

ttl = L"H_{eff}"*" resonances at "*L"N="*string(N)*", "*L"\Delta="*string(Δ)*", "*L"K="*string(K)*".\n Color bar in log scale. All levels considered."
heatmap(ϵ_1_array,ϵ_2_array, log.(matrix_color), xlab=L"\epsilon_1",ylab=L"\epsilon_2",title = ttl, guidefontsize=14)

# savefig("figs/important_figs/h_eff_resonances_heatmap_logscale.png")
