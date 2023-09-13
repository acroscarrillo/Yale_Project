using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs

df = DataFrame(CSV.File("data/h_eff_exact_resonances.csv"))

ttl = "Resonances at "*L"N="*string(N)*L",\Delta="*string(Δ)*L",K="*string(K)*",tol="*string(cross_tol)
xlb = L"\epsilon_1"
ylb = L"\epsilon_2"

# df_formatted = filter(row -> row.ΔE_n > 550 &&  row.ΔE_n < 650, df)  # match paper limits

# @df df_formatted plot(:ϵ_1, :ΔE_n, markersize = .5, xlab=L"\epsilon_1",ylab=L"E-E_0",seriestype=:scatter, legend=false,title=L"N="*string(N)*L",\Delta="*string(Δ)*L",K="*string(K)*L",\epsilon_2="*string(ϵ_2))


@df df plot(:ϵ_1, :ϵ_2, group=:cross_n,markersize = 2, xlab = xlb,ylab = ylb,seriestype = :scatter, title=ttl,alpha=1,guidefontsize=14,markerstrokewidth=0,legendtitle = "Levels", legendfontsize=10)


savefig("figs/h_eff_exact_resonances.png")