using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs

df = DataFrame(CSV.File("data/floquet_crossings_epsilon_1.csv"))
ttl = L"N="*string(N)*", "*L"\epsilon_2="*string(ϵ_2)*", "*L"\Delta="*string(Δ)

df_formatted = filter(row -> row.ΔE_n < 1000, df)  # match paper limits

df_floquet = filter(row -> row.floquet == 1, df_formatted)  # match paper limits
df_h_eff = filter(row -> row.floquet == 0, df_formatted) 
df_floquet.ϵ_1 .= (1/2.8).*(df_floquet.ϵ_1)

df_formatted = vcat(df_floquet,df_h_eff)


@df df_formatted plot(:ϵ_1, :ΔE_n,group=:floquet, markersize = 1, xlab=L"\epsilon_1/K",ylab=L"(E-E_0)/K,(\epsilon_n-\epsilon_0)/K",seriestype=:scatter, legend=true,title=ttl,markerstrokewidth=0,legendtitle="Floquet?",foreground_color_legend = nothing,background_color_legend=nothing)

@df df plot(:ϵ_1, :ΔE_n, group=:floquet,markersize = .5, xlab=L"\epsilon_1",ylab=L"\epsilon_n-\epsilon_0",seriestype=:scatter, legend=true,legendtitle=L"\epsilon_2",title=ttl,markerstrokewidth=0)

savefig("figs/floquet_crossing.png")