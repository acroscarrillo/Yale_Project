using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs

df = DataFrame(CSV.File("data/h_eff_repul.csv"))

df_formatted = filter(row -> row.ΔE_n > 550 &&  row.ΔE_n < 650, df)  # match paper limits

@df df_formatted plot(:ϵ_1, :ΔE_n, markersize = .5, xlab=L"\epsilon_1",ylab=L"E-E_0",seriestype=:scatter, legend=false,title=L"N="*string(N)*L",\Delta="*string(Δ)*L",K="*string(K)*L",\epsilon_2="*string(ϵ_2))

@df df plot(:ϵ_1, :ΔE_n, markersize = .5, xlab=L"\epsilon_1",ylab=L"E-E_0",seriestype=:scatter, legend=false,title=L"N="*string(N)*L",\Delta="*string(Δ)*L",K="*string(K)*L",\epsilon_2="*string(ϵ_2))

savefig("figs/h_eff_kissing.png")