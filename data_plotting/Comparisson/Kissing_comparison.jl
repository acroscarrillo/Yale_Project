using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs

df_comparison = DataFrame(CSV.File("data/kissing_comparison.csv"))

xlbl = L"ϵ_2/K"
ylbl = L"\Delta E_n/K,\Delta \epsilon _n/K" 
ttl = "g_n="*string(g_n)

#import comparison data and plot it
df_formatted = filter(row -> row.ΔE_n < 300, df_comparison)  # match paper limits
@df df_formatted plot(:ϵ_2, :ΔE_n, group=:Floquet, markersize = 1, xlab=xlbl, ylab=ylbl,seriestype=:scatter, legend=true,legendtitle="Floquet?",markerstrokewidth=0,title=ttl)

savefig("figs/kissing_comparison.png")