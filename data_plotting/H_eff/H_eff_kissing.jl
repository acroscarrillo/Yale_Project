using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs

df = DataFrame(CSV.File("data/h_eff_kissing.csv"))

df_formatted = filter(row -> row.ΔE_n > 0 &&  row.ΔE_n < 30*pi, df)  # match paper limits

@df df_formatted plot(:ξ, :ΔE_n, markersize = .5, xlab=L"\xi",ylab=L"E-E_0",seriestype=:scatter, legend=false)

savefig("figs/h_eff_kissing.png")