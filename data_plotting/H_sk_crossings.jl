include("../src/src.jl")

using StatsPlots
using LaTeXStrings

df = DataFrame(CSV.File("data/H_sk_spectrum.csv"))

df_formatted = filter(row -> row.E_n > 0 &&  row.E_n < 1400, df)

@df df_formatted plot(:ϵ_2_by_K, :E_n, markersize = 1, xlab=L"\epsilon_2 / K",ylab=L"E",seriestype=:scatter)


# @df df plot(:ϵ_2_by_K, -:E_n, markersize = 1, xlab=L"\epsilon_2 / K",ylab=L"E",seriestype=:scatter)




# savefig("figs/MI_vs_t/MI_dynamics_p_0_0.png")