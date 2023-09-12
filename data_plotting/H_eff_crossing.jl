include("../src/src.jl")

using StatsPlots
using LaTeXStrings

df = DataFrame(CSV.File("data/h_eff_crossing.csv"))

@df df contourf(:ϵ_1, :ϵ_2, :Δnn, markersize = .5, xlab=L"\epsilon_1",ylab=L"\epsilon_2",seriestype=:scatter, legend=false,title=L"N="*string(N)*L",K="*string(K)*L",\Delta="*string(Δ))

# @df df plot(:ϵ_2/K, :Δϵ_n, markersize = .5, xlab=L"\epsilon_2/K",ylab=L"\tilde{\epsilon} ",seriestype=:scatter, legend=false,title=L"N="*string(N)*L",\omega_0="*string(ω_0)*L",g_n="*string(g_n)*L",\omega_d="*string(ω_d),guidefontsize=15)

savefig("figs/h_eff_crossing.png")