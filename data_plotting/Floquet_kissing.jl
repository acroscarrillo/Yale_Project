include("../src/src.jl")

using StatsPlots
using LaTeXStrings

df = DataFrame(CSV.File("data/floquet_kissing.csv"))

@df df_floquet plot(:Ω_d, :Δϵ_n, markersize = .5, xlab=L"\Omega_d",ylab=L"\epsilon-\epsilon_0",seriestype=:scatter, legend=false,title=L"N="*string(N)*L",\omega_0="*string(ω_0)*L",g_n="*string(g_n)*L",\omega_d="*string(ω_d))

savefig("figs/floquet_kissing.png")