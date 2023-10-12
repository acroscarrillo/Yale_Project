using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs

df = DataFrame(CSV.File("data/comparison_crossings_epsilon_1.csv"))
ttl = L"N="*string(N)*", "*L"\epsilon_2/K="*string(round(ϵ_2_max,sigdigits=3))*", "*L"\Delta="*string(Δ)

df_formatted = filter(row -> row.ΔE_n < 400, df)  # match paper limits

# df_floquet = filter(row -> row.floquet == 1, df_formatted)  # match paper limits
# df_h_eff = filter(row -> row.floquet == 0, df_formatted) 
# df_floquet.ΔE_n .= (2).*(df_floquet.ΔE_n)
# df_floquet.ϵ_1 .= (2).*(df_floquet.ϵ_1)
# df_floquet.ϵ_2 .= (2).*(df_floquet.ϵ_2)
# df_formatted = vcat(df_floquet,df_h_eff)

@df df_formatted plot(:ϵ_1, :ΔE_n,group=:floquet, markersize = 1, xlab=L"\epsilon_1/K",ylab=L"(E-E_0)/K,(\epsilon_n-\epsilon_0)/K",seriestype=:scatter, legend=true,title=ttl,markerstrokewidth=0,legendtitle="Floquet?",foreground_color_legend = nothing,background_color_legend=nothing)

@df df plot(:ϵ_1, :ΔE_n,group=:floquet, markersize = 1, xlab=L"\epsilon_1/K",ylab=L"(E-E_0)/K,(\epsilon_n-\epsilon_0)/K",seriestype=:scatter, legend=true,title=ttl,markerstrokewidth=0,legendtitle="Floquet?",foreground_color_legend = nothing,background_color_legend=nothing)

savefig("figs/floquet_crossing.png")