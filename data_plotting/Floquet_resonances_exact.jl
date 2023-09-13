using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs

df = DataFrame(CSV.File("data/floquet_resonances.csv"))

ttl = "Resonances at "*L"N="*string(N)*L",\omega_0="*string(ω_0)*L",\omega_1="*string(ω_1)*L",\omega_2="*string(ω_2)*L",K="*string(K)*",tol="*string(cross_tol)

xlb = L"\epsilon_1"
ylb = L"\epsilon_2"


# @df df plot(:ϵ_1, :ϵ_2, group=:cross_n,markersize = 2, xlab = xlb,ylab = ylb,seriestype = :scatter, title=ttl,alpha=1,guidefontsize=14,markerstrokewidth=0,legendtitle = "Levels", legendfontsize=10)


@df df plot(:ϵ_1, :ϵ_2,markersize = 2, xlab = xlb,ylab = ylb,seriestype = :scatter, title=ttl,alpha=1,guidefontsize=14,markerstrokewidth=0, legend=false)


savefig("figs/h_eff_exact_resonances.png")