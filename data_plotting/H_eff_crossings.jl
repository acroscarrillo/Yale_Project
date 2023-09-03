include("../src/src.jl")

using StatsPlots
using LaTeXStrings

df = DataFrame(CSV.File("data/H_eff_spectrum.csv"))

df_formatted = filter(row -> row.ϵ_1 == 0 && row.K == -2.5, df)

@df df_formatted plot(:ϵ_2, :E_n, markersize = 1, xlab=L"\epsilon_2",ylab=L"E",seriestype=:scatter)





savefig("figs/MI_vs_t/MI_dynamics_p_0_0.png")