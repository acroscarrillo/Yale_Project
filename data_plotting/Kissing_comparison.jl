using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs

xlbl = L"ϵ_2/K"
ylbl = L"\Delta E_n/K,\Delta \epsilon _n/K" 

#import H_eff data and plot it
df_h_eff = DataFrame(CSV.File("data/h_eff_kissing.csv"))
df_formatted = filter(row -> row.ΔE_n < 400, df_h_eff)  # match paper limits
@df df_formatted plot(:ϵ_2, :ΔE_n, markersize = .5, xlab=xlbl, ylab=ylbl,seriestype=:scatter, legend=false)

#import Floquet data and plot it
df_floquet = DataFrame(CSV.File("data/floquet_kissing.csv"))
df_formatted = filter(row -> row.ΔE_n <400 && row.ΔE_n >0, df_floquet)  # match paper limits
@df df_formatted plot(:ϵ_2, :ΔE_n, markersize = .5, xlab=xlbl, ylab=ylbl,seriestype=:scatter, legend=false)
@df df_floquet plot(:ϵ_2, :ΔE_n, markersize = .5, xlab=xlbl, ylab=ylbl,seriestype=:scatter, legend=false)

#import comparison data and plot it
ttl = "Max and Rodri's ± issue.\\ MINUS. g_n="*string(g_n)
df_comparison = DataFrame(CSV.File("data/kissing_comparison.csv"))
df_formatted = filter(row -> row.ΔE_n < 300, df_comparison)  # match paper limits
@df df_formatted plot(:ϵ_2, :ΔE_n, group=:Floquet, markersize = 1, xlab=xlbl, ylab=ylbl,seriestype=:scatter, legend=true,legendtitle="Floquet?",markerstrokewidth=0,title=ttl)

@df df_comparison plot(:ϵ_2, :ΔE_n, group=:Floquet, markersize = 1, xlab=xlbl, ylab=ylbl,seriestype=:scatter, legend=true,legendtitle="Floquet?",markerstrokewidth=0)


savefig("figs/kissing_comparison.png")