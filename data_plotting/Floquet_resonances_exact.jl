using StatsPlots # plotting for idiots
using LaTeXStrings # latex support for sexy figs

# df = DataFrame(CSV.File("data/floquet_resonances.csv"))
df_floquet = DataFrame(CSV.File("data/data_to_gitignore/floquet_resonances_full_backup.csv"))
df_formatted = filter(row -> row.Δnn <= 0.000005, df_floquet)  # discard large crossings


ttl = "Resonances at "*L"N="*string(N)*L",\omega_0="*string(ω_0)*L",K="*string(K)#*",tol="*string(cross_tol)

ttl = ""

xlb = L"\epsilon_1"
ylb = L"\epsilon_2"

df_formatted = filter(row ->  row.n_cross < 5, df_formatted)  # match paper limits


@df df_formatted plot(:ϵ_1, :ϵ_2, group=:n_cross, markersize = 1, xlab = xlb,ylab = ylb,seriestype = :scatter,alpha=1,guidefontsize=14,markerstrokewidth=0, legend=true,legendtitle = "Levels", title = ttl,foreground_color_legend = nothing,background_color_legend=nothing,titlefontsize=10,legendpos=:right)

# @df df plot(:ϵ_1, :ϵ_2, group=:cross_n,markersize = 2, xlab = xlb,ylab = ylb,seriestype = :scatter, title=ttl,alpha=1,guidefontsize=14,markerstrokewidth=0,legendtitle = "Levels", legendfontsize=10)

df_formatted = filter(row ->  row.ϵ_2 < 0.0005, df)  # match paper limits

@df df_formatted plot(:ϵ_1, :ϵ_2, group=:n_cross, markersize = 1, xlab = xlb,ylab = ylb,seriestype = :scatter,alpha=1,guidefontsize=14,markerstrokewidth=0, legend=false,legendtitle = "Levels", title = ttl,foreground_color_legend = nothing,background_color_legend=nothing)



df_formatted = filter(row -> row.ϵ_2 == 4.500025267658904e-6, df)  # discard large crossings
@df df plot(:ϵ_1, :ϵ_2, group=:n_cross, markersize = 1, xlab = xlb,ylab = ylb,seriestype = :scatter,alpha=1,guidefontsize=14,markerstrokewidth=0, legend=true,legendtitle = "Levels", title = ttl,foreground_color_legend = nothing,background_color_legend=nothing)



# savefig("figs/important_figs/Floquet/floquet_resonances_heatmap_lin.png")